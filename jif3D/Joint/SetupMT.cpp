//============================================================================
// Name        : SetupMT.cpp
// Author      : Mar 1, 2010
// Version     : 
// Copyright   : 2010, mmoorkamp
//============================================================================

#include "SetupMT.h"
#include "../MT/ReadWriteImpedances.h"
#include "../MT/MTData.h"
#include "../MT/TipperData.h"
#include "../MT/X3DFieldCalculator.h"
#include "../Global/FileUtil.h"
#include "../Global/Noise.h"
#include "../Global/ReadWriteSparseMatrix.h"
#include "../Inversion/ModelTransforms.h"
#include "../Regularization/MinDiffRegularization.h"
#include "../ModelBase/ReadAnyModel.h"

#include <boost/make_shared.hpp>
#include <algorithm>
#include <string>

namespace jif3D
  {

    const std::string GridName = "Conductivity";
    const std::string DistName = "C";

    SetupMT::SetupMT() :
        GeneralDataSetup(GridName), relerr(0.02), tiprelerr(0.02), tipperr(0.0), DistCorr(
            0.0), tiplambda(0.0)
      {
      }

    SetupMT::~SetupMT()
      {
      }

    po::options_description SetupMT::SetupOptions()
      {
        //set the program option description for the MT part of the inversion
        po::options_description desc("MT options");
        desc.add_options()("mincond", po::value(&mincond)->default_value(1e-6))("maxcond",
            po::value(&maxcond)->default_value(5))("mtrelerr",
            po::value(&relerr)->default_value(0.02), "The relative error for the MT data")(
            "mtfine", po::value(&FineModelName),
            "The name for the model with the MT forward geometry")("mtinvcovar",
            po::value < std::string > (&MTInvCovarName),
            "Inverse covariance matrix to use in MT misfit calculation.")("inderrors",
            "Use the individual errors for each element instead of the same for all elements")(
            "x3dname", po::value < std::string > (&X3DName)->default_value("x3d"),
            "The name of the executable for x3d")("opt",
            "Use opt for Green's function calculation in x3d.")("distcorr",
            po::value(&DistCorr)->default_value(0),
            "Correct for distortion within inversion, value is regularization factor")(
            "tiprelerr", po::value(&tiprelerr)->default_value(0.02),
            "The relative error for the Tipper data")("tiperr",
            po::value(&tipperr)->default_value(0.01),
            "The absolute minimum error for the Tipper data")("cond_covmod",
            po::value<std::string>(),
            "A file containing the model covariance for conductivity");
        return desc;
      }

    bool SetupMT::SetupObjective(const po::variables_map &vm,
        jif3D::JointObjective &Objective, jif3D::ThreeDModelBase &InversionMesh,
        jif3D::rvec &CovModVec, std::vector<size_t> &startindices,
        std::vector<std::string> &SegmentNames, std::vector<parametertype> &SegmentTypes,
        boost::filesystem::path TempDir)
      {
        const size_t ngrid = InversionMesh.GetNModelElements();

        CovModVec.resize(ngrid);
        std::fill(CovModVec.begin(), CovModVec.end(), 1e-10);

        size_t startgrid = startindices.back();
        size_t endgrid = startgrid + ngrid;
        //first we ask the user a few questions
        //these are all values that are likely to change from run to run
        // so we do not want them as command line arguments or in the
        //configuration file
        std::cout << "MT Lambda: ";
        std::cin >> mtlambda;
        //if lambda is negative we assume that no MT inversion is wanted
        //and there is nothing to do here
        if (mtlambda > 0.0)
          {
            jif3D::MTData MTData;

            if (mtlambda > JointObjective::MinWeight)
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

                if (extension.compare(".dat") == 0)
                  {
                    MTData.ReadModEM(mtdatafilename);
                  }
                else
                  {
                    MTData.ReadNetCDF(mtdatafilename);
                  }
                MTData.WriteMeasurementPoints(mtdatafilename + ".mt_sites.vtk");
              }
            const size_t nmtsites = MTData.GetMeasPosX().size();

            const size_t ndata = MTData.GetData().size();
            std::string mtmodelfilename = jif3D::AskFilename("MT Model Filename: ");
            //read in the model and check whether the geometry matches the one
            //of the tomography starting model
            MTModel.ReadNetCDF(mtmodelfilename);
            if (mtlambda > JointObjective::MinWeight)
              {
                std::fill(CovModVec.begin(), CovModVec.end(), 1.0);

                //if x3d sees that a layer is completely homogeneous with the background
                //it optimizes this layer away which messes up our gradient calculation
                // as it does not output the field values for this layer
                //so we make the starting model slightly inhomogeneous to ensure
                //that this never happens
                for (size_t i = 0; i < MTModel.GetConductivities().shape()[2]; ++i)
                  {
                    MTModel.SetConductivities()[0][0][i] *= (1 + 0.0001 * (i + 1));
                  }

                bool WantDistCorr = (DistCorr > 0.0);
                //setup the objective function for the MT data
                boost::shared_ptr<jif3D::X3DFieldCalculator> FC = boost::make_shared<
                    jif3D::X3DFieldCalculator>(TempDir, X3DName);
                jif3D::X3DMTCalculator Calculator(TempDir, X3DName, WantDistCorr, true,
                    FC);

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
                        MinErr = jif3D::ConstructError(MTData.GetData(),
                            MTData.GetErrors(), relerr);
                      }
                    else
                      {
                        //use a relative value for the Berdichevskyi invariant at each period/site
                        MinErr = jif3D::ConstructMTError(MTData.GetData(), relerr);
                      }
                    //the error used in the inversion is the maximum of the error floor
                    //and the actual data error.
                    std::transform(MTData.GetErrors().begin(), MTData.GetErrors().end(),
                        MinErr.begin(), MinErr.begin(), [](double a, double b)
                          { return std::max(a,b);});
                    MTObjective->SetDataError(MinErr);
                  }


                boost::shared_ptr<jif3D::ChainedTransform> ConductivityTransform(
                    new jif3D::ChainedTransform);
                jif3D::rvec RefModel(ngrid);
                std::fill(RefModel.begin(), RefModel.end(), 1.0);
                ConductivityTransform->AppendTransform(
                    boost::shared_ptr<jif3D::GeneralModelTransform>(
                        new jif3D::TanhTransform(std::log(mincond), std::log(maxcond))));
                ConductivityTransform->AppendTransform(
                    boost::shared_ptr<jif3D::GeneralModelTransform>(
                        new jif3D::LogTransform(RefModel)));


                //add the MT part to the JointObjective that will be used
                //for the inversion
                Transform = boost::make_shared<jif3D::MultiSectionTransform>(2 * ngrid,
                    startgrid, endgrid, ConductivityTransform);

                //output some information to the screen
                //to signal that we added the MT data
                std::cout << "MT ndata: " << ndata << std::endl;
                std::cout << "MT lambda: " << mtlambda << std::endl;

                if (vm.count("cond_covmod"))
                  {
                    boost::shared_ptr<jif3D::ThreeDModelBase> CovModel =
                        jif3D::ReadAnyModel(vm["cond_covmod"].as<std::string>());

                    const size_t ncovmod = CovModel->GetData().num_elements();
                    if (ncovmod != ngrid)
                      {
                        throw jif3D::FatalException(
                            " Conductivity covariance does not have the same number of parameters "
                                + std::to_string(ncovmod) + " as inversion model "
                                + std::to_string(ngrid), __FILE__, __LINE__);
                      }
                    else
                      {
                        std::cout << " Setting covariance for conductivity from file. "
                            << std::endl;
                      }
                    std::copy(CovModel->GetData().origin(),
                        CovModel->GetData().origin() + ncovmod, CovModVec.begin());

                  }
                //if we want to correct for distortion within the inversion
                boost::shared_ptr<jif3D::GeneralModelTransform> Copier(
                    new jif3D::ModelCopyTransform);
                if (WantDistCorr)
                  {

                    jif3D::rvec CRef(nmtsites * 4);
                    for (size_t i = 0; i < nmtsites; ++i)
                      {
                        CRef(i * 4) = 1.0;
                        CRef(i * 4 + 1) = 0.0;
                        CRef(i * 4 + 2) = 0.0;
                        CRef(i * 4 + 3) = 1.0;
                      }
                    size_t startdist = endgrid;
                    size_t enddist = startdist + CRef.size();

                    //this transformation only becomes active if we use distortion correction with MT data
                    //in this case we add extra inversion parameters beyond the current ones
                    //so we set it up here that we can access it later, but with a parameter setting that only
                    //works if we actually have distortion correction, so we have to be careful later
                    auto DistRegTrans = boost::make_shared<jif3D::MultiSectionTransform>(
                        enddist, startdist, enddist, Copier);

                    std::vector<double> C = MTData.GetDistortion();

                    jif3D::X3DModel DistModel;
                    DistModel.SetMeshSize(nmtsites * 4, 1, 1);
                    boost::shared_ptr<jif3D::RegularizationFunction> DistReg(
                        new jif3D::MinDiffRegularization(DistModel));
                    DistReg->SetReferenceModel(CRef);

                    Objective.AddObjective(DistReg, DistRegTrans, DistCorr, "DistReg",
                        jif3D::JointObjective::regularization);
                    startindices.push_back(enddist);
                    SegmentNames.push_back(DistName);
                    SegmentTypes.push_back(GeneralDataSetup::other);

                    Transform->AddSection(ngrid, ngrid + CRef.size(), Copier);
                    StartingParameters.resize(ngrid + CRef.size());
                    std::copy(MTModel.GetConductivities().origin(),
                        MTModel.GetConductivities().origin() + ngrid,
                        StartingParameters.begin());
                    std::copy(C.begin(), C.end(), StartingParameters.begin() + ngrid);

                  }
                else
                  {
                    StartingParameters.resize(ngrid);
                    std::copy(MTModel.GetConductivities().origin(),
                        MTModel.GetConductivities().origin() + ngrid,
                        StartingParameters.begin());
                  }

                Objective.AddObjective(MTObjective, Transform, mtlambda, "MT",
                    JointObjective::datafit);
                //ask for Tipper, currently it only works in conjunction with MT
                //otherwise we would need to check for x3d and conductivity models

                std::cout << "Tipper Lambda: ";
                std::cin >> tiplambda;
                if (tiplambda > 0)
                  {
                    jif3D::TipperData DataTip;
                    std::string tipdatafilename = jif3D::AskFilename(
                        "Tipper data filename: ");
                    std::string extension = jif3D::GetFileExtension(tipdatafilename);
                    DataTip.ReadNetCDF(tipdatafilename);

                    jif3D::X3DTipperCalculator TipCalc(TempDir, X3DName, true, FC);
                    TipObjective = boost::make_shared<
                        jif3D::ThreeDModelObjective<jif3D::X3DTipperCalculator>>(TipCalc);
                    TipObjective->SetCoarseModelGeometry(MTModel);
                    TipObjective->SetObservedData(DataTip);

                    size_t ntip = DataTip.GetData().size();
                    jif3D::rvec TE(ntip, 0.0);
                    std::vector<double> TError = DataTip.GetErrors();
                    for (size_t i = 0; i < ntip / 4; ++i)
                      {
                        double abs_re = std::sqrt(
                            std::pow(DataTip.GetData().at(4 * i), 2)
                                + std::pow(DataTip.GetData().at(4 * i + 2), 2));
                        double abs_im = std::sqrt(
                            std::pow(DataTip.GetData().at(4 * i + 1), 2)
                                + std::pow(DataTip.GetData().at(4 * i + 3), 2));
                        TE(4 * i) = std::max(
                            std::max(abs_re * tiprelerr, TError.at(4 * i)),
                            TError.at(4 * i + 2));
                        TE(4 * i + 2) = TE(4 * i);
                        TE(4 * i + 1) = std::max(
                            std::max(abs_im * tiprelerr, TError.at(4 * i)),
                            TError.at(4 * i + 2));
                        TE(4 * i + 3) = TE(4 * i + 1);
                      }
                    std::copy(TE.begin(), TE.end(), TError.begin());

                    std::transform(TError.begin(), TError.end(), TError.begin(),
                        [this](double a)
                          { return std::max(a,tipperr);});

                    TipObjective->SetDataError(TError);
                    Objective.AddObjective(TipObjective, Transform, tiplambda, "Tipper",
                        jif3D::JointObjective::datafit);

                  }
              }


            startindices.push_back(endgrid);
            SegmentNames.push_back(GridName);
            SegmentTypes.push_back(GeneralDataSetup::gridparameter);
          }
        //return true if we added an MT objective function
        return (mtlambda > JointObjective::MinWeight);
      }

    void SetupMT::IterationOutput(const std::string &filename,
        const jif3D::rvec &ModelVector)
      {
        if (mtlambda > JointObjective::MinWeight)
          {
            jif3D::rvec CondInvModel = Transform->GeneralizedToPhysical(ModelVector);
            const size_t ngrid = MTModel.GetNModelElements();
            std::copy(CondInvModel.begin(), CondInvModel.begin() + ngrid,
                MTModel.SetConductivities().origin());
            MTModel.WriteVTK(filename + ".mt.inv.vtk");
            MTModel.WriteNetCDF(filename + ".mt.inv.nc");
          }
      }

    void SetupMT::FinalOutput(const std::string &filename,
        const jif3D::rvec &FinalModelVector)
      {
        if (mtlambda > JointObjective::MinWeight)
          {
            std::cout << "Writing final conductivity models " << std::endl;
            const size_t ngrid = MTModel.GetNModelElements();
            jif3D::rvec CondInvModel = Transform->GeneralizedToPhysical(FinalModelVector);
            std::copy(CondInvModel.begin(), CondInvModel.begin() + ngrid,
                MTModel.SetConductivities().origin());

            auto MTData = MTObjective->GetObservedData();
            auto MTDataVec = MTObjective->GetSyntheticData();
            auto MTErrVec = MTObjective->GetDataError();
            if (MTDataVec.size() > 0)
              {
                MTData.SetDataAndErrors(
                    std::vector<double>(MTDataVec.begin(), MTDataVec.end()), MTErrVec);
                MTData.WriteNetCDF(filename + ".inv_imp.nc");
              }
            MTModel.WriteVTK(filename + ".mt.inv.vtk");
            MTModel.WriteNetCDF(filename + ".mt.inv.nc");
          }
      }
  }
