//============================================================================
// Name        : SetupMT.cpp
// Author      : Mar 1, 2010
// Version     : 
// Copyright   : 2010, mmoorkamp
//============================================================================

#include "../MT/ReadWriteImpedances.h"
#include "../MT/MTData.h"
#include "../Global/FileUtil.h"
#include "../Global/Noise.h"
#include "../Global/ReadWriteSparseMatrix.h"

#include "SetupMT.h"
#include <boost/make_shared.hpp>
#include <algorithm>
#include <string>

namespace jif3D
  {

    SetupMT::SetupMT() :
        relerr(0.02), DistCorr(0.0)
      {
      }

    SetupMT::~SetupMT()
      {
      }

    po::options_description SetupMT::SetupOptions()
      {
        //set the program option description for the MT part of the inversion
        po::options_description desc("MT options");
        desc.add_options()("mtrelerr", po::value(&relerr)->default_value(0.02),
            "The relative error for the MT data")("mtfine", po::value(&FineModelName),
            "The name for the model with the MT forward geometry")("mtinvcovar",
            po::value<std::string>(&MTInvCovarName),
            "Inverse covariance matrix to use in MT misfit calculation.")("inderrors",
            "Use the individual errors for each element instead of the same for all elements")(
            "x3dname", po::value<std::string>(&X3DName)->default_value("x3d"),
            "The name of the executable for x3d")("opt",
            "Use opt for Green's function calculation in x3d.")("distcorr",
            po::value(&DistCorr)->default_value(0),
            "Correct for distortion within inversion, value is regularization factor");
        return desc;
      }

    bool SetupMT::SetupObjective(const po::variables_map &vm,
        jif3D::JointObjective &Objective,
        boost::shared_ptr<jif3D::GeneralModelTransform> Transform, double xorigin,
        double yorigin, boost::filesystem::path TempDir)
      {
        //first we ask the user a few questions
        //these are all values that are likely to change from run to run
        // so we do not want them as command line arguments or in the
        //configuration file
        double mtlambda = 1.0;
        std::cout << "MT Lambda: ";
        std::cin >> mtlambda;
        //if lambda is negative we assume that no MT inversion is wanted
        //and there is nothing to do here
        if (mtlambda > 0.0)
          {
            if (!boost::filesystem::exists(X3DName))
              {
                std::cerr << X3DName << " is not accessible or  does not exist ! \n";
                return 500;
              }

            //for inversion we need some data, so we ask for the filename
            std::string mtdatafilename = jif3D::AskFilename("MT data filename: ");
            std::string extension = jif3D::GetFileExtension(mtdatafilename);
            //read in MT data, the position of the measurement sites, frequencies and impedances
            // we also try to read in the parameters of the distortion Matrix C
            //if these are not present they will be set to identity matrix for each site
            //in the forward calculation, otherwise the synthetic responses will be multiplied

            jif3D::MTData MTData;

            if (extension.compare(".dat") == 0)
              {
                MTData.ReadModEM(mtdatafilename);
              }
            else
              {
                MTData.ReadNetCDF(mtdatafilename);
              }
            const size_t ndata = MTData.GetData().size();
            std::string mtmodelfilename = jif3D::AskFilename("MT Model Filename: ");
            //read in the model and check whether the geometry matches the one
            //of the tomography starting model
            MTModel.ReadNetCDF(mtmodelfilename);

            //if x3d sees that a layer is completely homogeneous with the background
            //it optimizes this layer away which messes up our gradient calculation
            // as it does not output the field values for this layer
            //so we make the starting model slightly inhomogeneous to ensure
            //that this never happens
            for (size_t i = 0; i < MTModel.GetConductivities().shape()[2]; ++i)
              {
                MTModel.SetConductivities()[0][0][i] *= (1 + 0.0001 * (i + 1));
              }
            if (xorigin != 0.0 || yorigin != 0.0)
              {
                MTModel.SetOrigin(xorigin, yorigin, 0.0);
              }
            bool WantDistCorr = (DistCorr > 0.0);
            //setup the objective function for the MT data
            jif3D::X3DMTCalculator Calculator(TempDir, X3DName, WantDistCorr);

            if (vm.count("opt"))
              {
                std::cout << "Using Opt type Green's functions " << std::endl;
                Calculator.SetGreenType1(jif3D::GreenCalcType::opt);
                Calculator.SetGreenType4(jif3D::GreenCalcType::opt);
              }

            MTObjective = boost::make_shared<
                jif3D::ThreeDModelObjective<jif3D::X3DMTCalculator> >(Calculator);
            //if we specified the name for a refined model for forward calculations
            //we read in that model, set the measurement configuration for the observed
            //data and pass it to the objective function
            if (vm.count("mtfine"))
              {
                jif3D::X3DModel FineModel;
                FineModel.ReadNetCDF(FineModelName);
                MTObjective->SetFineModelGeometry(FineModel);
              }
            MTObjective->SetCoarseModelGeometry(MTModel);
            MTObjective->SetObservedData(MTData);
            if (vm.count("mtinvcovar"))
              {
                comp_mat InvCov(ndata, ndata);
                ReadSparseMatrixFromNetcdf(MTInvCovarName, InvCov, "InvCovariance");
                MTObjective->SetInvCovMat(InvCov);
              }
            else
              {
                //construct an error floor for each impedance element
                std::vector<double> MinErr(MTData.GetErrors().size());
                if (vm.count("inderrors"))
                  {
                    //use a relative value for each impedance element separately
                    MinErr = jif3D::ConstructError(MTData.GetData(), MTData.GetErrors(), relerr);
                  }
                else
                  {
                    //use a relative value for the Berdichevskyi invariant at each period/site
                    MinErr = jif3D::ConstructMTError(MTData.GetData(), relerr);
                  }
                //the error used in the inversion is the maximum of the error floor
                //and the actual data error.
                std::transform(MTData.GetErrors().begin(),MTData.GetErrors().end(), MinErr.begin(),
                    MinErr.begin(), [](double a, double b)
                      { return std::max(a,b);});
                MTObjective->SetDataError(MinErr);
              }
            //add the MT part to the JointObjective that will be used
            //for the inversion
            Objective.AddObjective(MTObjective, Transform, mtlambda, "MT",
                JointObjective::datafit);
            //output some information to the screen
            //to signal that we added the MT data
            std::cout << "MT ndata: " << ndata << std::endl;
            std::cout << "MT lambda: " << mtlambda << std::endl;
          }
        //return true if we added an MT objective function
        return (mtlambda > 0.0);
      }
  }
