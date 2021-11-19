//============================================================================
// Name        : SetupMagnetics.cpp
// Author      : Mar 1, 2010
// Version     : 
// Copyright   : 2010, mmoorkamp
//============================================================================

#include "SetupMagnetics.h"
#include "../Inversion/ModelTransforms.h"
#include "../Magnetics/ReadWriteMagneticData.h"
#include "../Magnetics/OMPMagneticImp.h"
#include "../Magnetics/MagneticTransforms.h"
#include "../Gravity/DepthWeighting.h"
#include "../Global/FileUtil.h"
#include "../Global/Noise.h"
#include <boost/make_shared.hpp>

namespace jif3D
  {

    SetupMagnetics::SetupMagnetics() :
        inclination(0.0), declination(0.0), fieldstrength(1.0), relerr(0.02), minerr(0.0)
      {

      }

    SetupMagnetics::~SetupMagnetics()
      {

      }

    po::options_description SetupMagnetics::SetupOptions()
      {
        po::options_description desc("Magnetics options");

        desc.add_options()("gpu", "Perform Magnetics calculation on GPU")("magrelerr",
            po::value(&relerr)->default_value(0.02),
            "The relative error for the Magnetics data")("magminerr",
            po::value(&minerr)->default_value(0.0),
            "The minimum absolute error for the Magnetics data")("inclination",
            po::value(&inclination)->default_value(0.0),
            "The inclination of the magnetic field in radian")("declination",
            po::value(&declination)->default_value(0.0),
            "The declination of the magnetic field in radian")("fieldstrength",
            po::value(&fieldstrength)->default_value(1.0),
            "The strength of the inducing magnetic field in nT");

        return desc;
      }

    bool SetupMagnetics::SetupObjective(const po::variables_map &vm,
        jif3D::JointObjective &Objective,
        boost::shared_ptr<jif3D::GeneralModelTransform> &Transform, double xorigin,
        double yorigin, boost::filesystem::path TempDir)
      {
        //if we want to use CUDA for forward modeling
        //we set a variable for easier access later and print
        //that information to the screen
        bool wantcuda = false;
        if (vm.count("gpu"))
          {
            std::cout << "Using GPU" << "\n";
            wantcuda = true;
          }

        jif3D::MagneticData MagData;
        double maglambda = 1.0;

        //we first ask for the weights for scalar and tensor Magnetics
        //as we might not have to do very much if the weights are zero
        std::cout << "Magnetics Lambda: ";
        std::cin >> maglambda;

        //if the weight is different from zero
        //we have to deal with magnetics
        if (maglambda > 0.0)
          {
            //if the weight is really small we use this as a "secret setting"
            //for constrained inversion, so we only read in data when the
            //weight is larger than 1e-32, otherwise we create a single point
            //dummy dataset, the joint objective class is set that
            //objective functions with such a small weight are not evaluated anyway
            if (maglambda > JointObjective::MinWeight)
              {
                //read in data file
                std::string magdatafilename = jif3D::AskFilename(
                    "Total field magnetic Data Filename: ");
                MagData.ReadNetCDF(magdatafilename);
                MagData.WriteVTK(magdatafilename);
              }
            else
              {
                //create dummy dataset to make constrained inversion more convenient
                //otherwise the user would have to do this or we would have to create extra settings
                MagData.AddMeasurementPoint(0.0, 0.0, 0.0);
                std::vector<double> dummy(1, 1.0);
                MagData.SetDataAndErrors(dummy, dummy);
              }
            std::string magmodelfilename = jif3D::AskFilename(
                "Magnetics Model Filename: ");
            Model.ReadNetCDF(magmodelfilename);
            if (xorigin != 0.0 || yorigin != 0.0)
              {
                Model.SetOrigin(xorigin, yorigin, 0.0);
              }

            //we want to set the path for temporary file storage
            //the factory function cannot perform this, so we
            //have to assemble the calculator object ourselves
            boost::shared_ptr<jif3D::ThreeDGravMagImplementation<jif3D::MagneticData> > Implementation;
            if (wantcuda)
              {
                throw jif3D::FatalException("No GPU support, yet !", __FILE__, __LINE__);
              }
            else
              {
                Implementation = boost::shared_ptr<
                    jif3D::ThreeDGravMagImplementation<jif3D::MagneticData> >(
                    new jif3D::OMPMagneticImp(inclination, declination, fieldstrength));
              }
#ifdef MAGDISK
            Calculator = boost::make_shared<MagCalculatorType>(Implementation, TempDir);
#else
            Calculator = boost::make_shared<MagCalculatorType>(Implementation);
#endif
            Calculator->SetDataTransform(
                boost::shared_ptr<jif3D::TotalFieldAnomaly>(
                    new jif3D::TotalFieldAnomaly(inclination, declination)));

            MagObjective =
                boost::shared_ptr<jif3D::ThreeDModelObjective<MagCalculatorType> >(
                    new jif3D::ThreeDModelObjective<MagCalculatorType>(*Calculator));
            MagObjective->SetObservedData(MagData);
            MagObjective->SetCoarseModelGeometry(Model);
            std::vector<double> Error(
                jif3D::ConstructError(MagData.GetData(), MagData.GetErrors(), relerr,
                    minerr));
            MagObjective->SetDataError(Error);

            Objective.AddObjective(MagObjective, Transform, maglambda, "Magnetics",
                JointObjective::datafit);
            std::cout << " Magnetics ndata: " << MagData.GetData().size() << std::endl;
            std::cout << " Magnetics lambda: " << maglambda << std::endl;
          }

        //indicate whether we added a Magnetics objective function
        //this way the caller can do additional consistency checks
        return (maglambda > JointObjective::MinWeight);
      }
  }
