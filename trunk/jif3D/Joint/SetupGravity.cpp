//============================================================================
// Name        : SetupGravity.cpp
// Author      : Mar 1, 2010
// Version     : 
// Copyright   : 2010, mmoorkamp
//============================================================================

#include "SetupGravity.h"
#include "../Inversion/ModelTransforms.h"
#include "../Gravity/ReadWriteGravityData.h"
#include "../Gravity/ThreeDGravityFactory.h"
#include "../Gravity/ScalarGravityData.h"
#include "../Gravity/TensorGravityData.h"
#include "../Gravity/ReadWriteGravityData.h"
#include "../Gravity/DepthWeighting.h"

#include "../ModelBase/VTKTools.h"
#include "../ModelBase/ThreeDModelBase.h"
#include "../ModelBase/ReadAnyModel.h"

#include "../Global/FileUtil.h"
#include "../Global/Noise.h"
#include "../Inversion/ModelTransforms.h"

namespace jif3D
  {

    const std::string Name = "Density";
    SetupGravity::SetupGravity() :
        GeneralDataSetup(Name), scalrelerr(0.02), ftgrelerr(0.02), scalminerr(0.0), ftgminerr(
            1e-9), HaveScal(false), HaveFTG(false)
      {

      }

    SetupGravity::~SetupGravity()
      {

      }

    po::options_description SetupGravity::SetupOptions()
      {
        po::options_description desc("Gravity options");

        desc.add_options()("mindens", po::value(&mindens)->default_value(-1000.0))(
            "maxdens", po::value(&maxdens)->default_value(1000.0))("gpu",
            "Perform gravity calculation on GPU")("scalrelerr",
            po::value(&scalrelerr)->default_value(0.01),
            "The relative error for the scalar gravity data")("scalminerr",
            po::value(&scalminerr)->default_value(1e-5),
            "The minimum absolute error for the scalar gravity data")("ftgrelerr",
            po::value(&ftgrelerr)->default_value(0.01),
            "The relative error for the FTG gravity data")("ftgminerr",
            po::value(&ftgminerr)->default_value(1e-9),
            "The minimum absolute error for the FTG gravity data")("dens_covmod",
            po::value<std::string>(),
            "A file containing the model covariance for densities")("gravdepth",
            "Counteract the decay in sensitivities of gravity data with depth");

        return desc;
      }

    bool SetupGravity::SetupObjective(const po::variables_map &vm,
        jif3D::JointObjective &Objective, jif3D::ThreeDModelBase &InversionMesh,
        jif3D::rvec &CovModVec, std::vector<size_t> &startindices,
        std::vector<std::string> &SegmentNames, std::vector<parametertype> &SegmentTypes,
        boost::filesystem::path TempDir)
      {
        //if we want to use CUDA for forward modeling
        //we set a variable for easier access later and print
        //that information to the screen
        const size_t ngrid = InversionMesh.GetNModelElements();
        bool wantcuda = false;
        if (vm.count("gpu"))
          {
            std::cout << "Using GPU" << "\n";
            wantcuda = true;
          }

        //we initially assume that we do not want to change the density model
        //but use it only as a constraint, if the weight is below the minimal threshold
        //and we read in data, we set it to one below, also it gets overwritten if we specify a covariance file
        CovModVec.resize(ngrid);
        std::fill(CovModVec.begin(), CovModVec.end(), 1e-10);

        jif3D::ScalarGravityData ScalGravData;
        jif3D::TensorGravityData FTGData;
        double scalgravlambda = 1.0;
        double ftglambda = 1.0;
        //we first ask for the weights for scalar and tensor gravity
        //as we might not have to do very much if the weights are zero
        std::cout << "Scalar Gravimetry Lambda: ";
        std::cin >> scalgravlambda;
        std::cout << "FTG Lambda: ";
        std::cin >> ftglambda;

        //if the weight is different from zero
        //we have to read in scalar gravity data
        if (scalgravlambda > JointObjective::MinWeight)
          {
            std::string scalgravdatafilename = jif3D::AskFilename(
                "Scalar Gravity Data Filename: ");

            ScalGravData.ReadNetCDF(scalgravdatafilename);
            ScalGravData.WriteVTK(scalgravdatafilename);
            HaveScal = true;
            std::fill(CovModVec.begin(), CovModVec.end(), 1.0);

          }

        //if the weight is different from zero
        //we have to read in ftg data
        if (ftglambda > JointObjective::MinWeight)
          {
            std::string ftgdatafilename = jif3D::AskFilename("FTG Data Filename: ");
            FTGData.ReadNetCDF(ftgdatafilename);
            HaveFTG = true;
            std::fill(CovModVec.begin(), CovModVec.end(), 1.0);
          }
        //if the inversion includes any type of gravity data
        //we need the model geometry
        if (scalgravlambda > 0.0 || ftglambda > 0.0)
          {
            std::string gravmodelfilename = jif3D::AskFilename(
                "Gravity Model Filename: ");
            ScalGravModel.ReadNetCDF(gravmodelfilename);

            if (ngrid != ScalGravModel.GetDensities().num_elements())
              {
                std::string err = "Gravity model does not match grid size ! "
                    + std::to_string(ScalGravModel.GetDensities().num_elements()) + " "
                    + std::to_string(ngrid);
                throw jif3D::FatalException(err, __FILE__, __LINE__);
              }
            StartingParameters.resize(ngrid);
            std::copy(ScalGravModel.GetDensities().origin(),
                ScalGravModel.GetDensities().origin() + ngrid,
                StartingParameters.begin());

            if (vm.count("dens_covmod"))
              {
                boost::shared_ptr<jif3D::ThreeDModelBase> CovModel = jif3D::ReadAnyModel(
                    vm["dens_covmod"].as<std::string>());

                const size_t ncovmod = CovModel->GetData().num_elements();
                if (ncovmod != ngrid)
                  {
                    std::cerr
                        << " Density covariance does not have the same number of parameters "
                        << ncovmod << " as inversion model " << ngrid << std::endl;
                    return 100;
                  }
                else
                  {
                    std::cout << " Setting covariance for density from file. "
                        << std::endl;
                  }
                std::copy(CovModel->GetData().origin(),
                    CovModel->GetData().origin() + ncovmod, CovModVec.begin());
              }

            if (vm.count("gravdepth"))
              {

                boost::shared_ptr<
                    jif3D::ThreeDGravMagImplementation<jif3D::ScalarGravityData> > Implementation =
                    boost::make_shared<jif3D::ScalarOMPGravityImp>();
                jif3D::FullSensitivityGravMagCalculator<jif3D::ScalarGravityData> FullCalc(
                    Implementation);

                std::cout << "Calculating depth weighting." << std::endl;
                //now we perform the depth weighting for the sensitivities
                jif3D::rvec SensProfile, WeightVector;
                //we find a measurement site close to the centre of the model and extract the
                //sensitivity variation with depth
                jif3D::rmat Sens;
                jif3D::CalculateMiddleSens(ScalGravModel, FullCalc, SensProfile);

                double DepthExponent = -2.0;
                //we fit a curve of the form 1/(z+z0)^n to the extracted sensitivities
                double z0 = FitZ0(SensProfile, ScalGravModel.GetZCellSizes(),
                    jif3D::WeightingTerm(DepthExponent));
                std::cout << "Estimated z0: " << z0 << std::endl;
                const size_t zsize = ScalGravModel.GetModelShape()[2];
                //calculate the depth scaling
                jif3D::ConstructDepthWeighting(ScalGravModel.GetZCellSizes(), z0,
                    WeightVector, jif3D::WeightingTerm(DepthExponent));
                for (size_t i = 0; i < ngrid; ++i)
                  {
                    CovModVec(i) *= WeightVector(i % zsize);
                  }
                std::ofstream proffile("profile.out");
                std::copy(SensProfile.begin(), SensProfile.end(),
                    std::ostream_iterator<double>(proffile, "\n"));
                std::ofstream covfile("cov.out");
                std::copy(CovModVec.begin(), CovModVec.end(),
                    std::ostream_iterator<double>(covfile, "\n"));
                std::ofstream weightfile("weights.out");
                std::copy(WeightVector.begin(), WeightVector.end(),
                    std::ostream_iterator<double>(weightfile, "\n"));
                jif3D::ThreeDGravityModel DepthModel(ScalGravModel);
                std::copy(CovModVec.begin(), CovModVec.begin() + ngrid,
                    DepthModel.SetDensities().origin());
                DepthModel.WriteNetCDF("depth_grav_cov.nc");
                DepthModel.WriteVTK("depth_grav_cov.vtk");
              }

            boost::shared_ptr<jif3D::GeneralModelTransform> DensTrans =
                boost::make_shared<jif3D::TanhTransform>(mindens, maxdens);
            size_t start = startindices.back();
            size_t end = start + ngrid;
            Transform = boost::make_shared<jif3D::MultiSectionTransform>(2 * ngrid, start,
                end, DensTrans);
            startindices.push_back(end);
            SegmentNames.push_back(Name);
            SegmentTypes.push_back(GeneralDataSetup::gridparameter);
          }
        //now we setup the objective functions for each type of gravity data
        //the steps are the same for both, create new objective function,
        //set the observed data, model geometry and data error
        //finally output some basic information about the data to the screen
        //to signal the user that something happened.
        if (scalgravlambda > JointObjective::MinWeight)
          {
            //we want to set the path for temporary file storage
            //the factory function cannot perform this, so we
            //have to assemble the calculator object ourselves
            boost::shared_ptr<jif3D::ThreeDGravMagImplementation<jif3D::ScalarGravityData> > Implementation;
            if (wantcuda)
              {
#ifdef HAVEGPU
                Implementation = boost::shared_ptr<jif3D::ThreeDGravMagImplementation<jif3D::ScalarGravityData> >(
                    new jif3D::ScalarCudaGravityImp);
#else
                throw jif3D::FatalException(
                    "Code has been compiled without GPU support !", __FILE__, __LINE__);
#endif
              }
            else
              {
                Implementation = boost::shared_ptr<
                    jif3D::ThreeDGravMagImplementation<jif3D::ScalarGravityData> >(
                    new jif3D::ScalarOMPGravityImp);
              }
#ifdef GRAVDISK
            ScalarCalculatorType ScalarCalculator(Implementation, TempDir);
            std::cout << "Scalar gravity will take " << ngrid * ScalGravData.GetData().size() * 8 / 1e9 << " GB disk space " << std::endl;
#else
            ScalarCalculatorType ScalarCalculator(Implementation);
#endif
            ScalGravObjective = boost::make_shared<
                jif3D::ThreeDModelObjective<ScalarCalculatorType> >(ScalarCalculator);
            ScalGravObjective->SetObservedData(ScalGravData);
            ScalGravObjective->SetCoarseModelGeometry(ScalGravModel);
            ScalGravObjective->SetDataError(
                jif3D::ConstructError(ScalGravData.GetData(), ScalGravData.GetErrors(),
                    scalrelerr, scalminerr));

            Objective.AddObjective(ScalGravObjective, Transform, scalgravlambda,
                "ScalGrav", JointObjective::datafit);
            std::cout << "Scalar Gravity ndata: " << ScalGravData.GetData().size()
                << std::endl;
            std::cout << "Scalar Gravity lambda: " << scalgravlambda << std::endl;
          }
        if (ftglambda > JointObjective::MinWeight)
          {
            //we want to set the path for temporary file storage
            //the factory function cannot perform this, so we
            //have to assemble the calculator object ourselves
            boost::shared_ptr<jif3D::ThreeDGravMagImplementation<jif3D::TensorGravityData> > Implementation;
            if (wantcuda)
              {
#ifdef HAVEGPU
                Implementation = boost::shared_ptr<jif3D::ThreeDGravMagImplementation<jif3D::TensorGravityData> >(
                    new jif3D::TensorCudaGravityImp);
#else
                throw jif3D::FatalException(
                    "Code has been compiled without GPU support !", __FILE__, __LINE__);
#endif
              }
            else
              {
                Implementation = boost::shared_ptr<
                    jif3D::ThreeDGravMagImplementation<jif3D::TensorGravityData> >(
                    new jif3D::TensorOMPGravityImp);
              }
#ifdef GRAVDISK
            TensorCalculatorType TensorCalculator(Implementation, TempDir);
            std::cout << "Tensor gravity will take " << ngrid * ScalGravData.GetData().size() * 8 * 9/ 1e9 << " GB disk space " << std::endl;
#else
            TensorCalculatorType TensorCalculator(Implementation);
#endif
            FTGObjective = boost::make_shared<
                jif3D::ThreeDModelObjective<TensorCalculatorType> >(TensorCalculator);
            FTGObjective->SetObservedData(FTGData);
            FTGObjective->SetCoarseModelGeometry(ScalGravModel);
            FTGObjective->SetDataError(
                jif3D::ConstructError(FTGData.GetData(), FTGData.GetErrors(), ftgrelerr,
                    ftgminerr));

            Objective.AddObjective(FTGObjective, Transform, ftglambda, "FTG",
                JointObjective::datafit);
            std::cout << "FTG ndata: " << FTGData.GetData().size() << std::endl;
            std::cout << "FTG lambda: " << ftglambda << std::endl;
          }
        //indicate whether we added a gravity objective function
        //this way the caller can do additional consistency checks
        return (ftglambda > JointObjective::MinWeight)
            || (scalgravlambda > JointObjective::MinWeight);
      }

    void SetupGravity::IterationOutput(const std::string &filename,
        const jif3D::rvec &ModelVector)
      {
        if (HaveScal || HaveFTG)
          {
            const size_t ngrid = ScalGravModel.GetNModelElements();

            jif3D::rvec DensVec = Transform->GeneralizedToPhysical(ModelVector);
            std::copy(DensVec.begin(), DensVec.begin() + ngrid,
                ScalGravModel.SetDensities().origin());
            ScalGravModel.WriteVTK(filename + ".grav.inv.vtk");
            ScalGravModel.WriteNetCDF(filename + ".grav.inv.nc");

          }
      }

    void SetupGravity::FinalOutput(const std::string &filename ,const jif3D::rvec &FinalModelVector)
      {

        const size_t ngrid = ScalGravModel.GetNModelElements();

        if (HaveScal || HaveFTG)
          {
            std::cout << "Writing final density models " << std::endl;
            jif3D::rvec DensVec = Transform->GeneralizedToPhysical(FinalModelVector);
            std::copy(DensVec.begin(), DensVec.begin() + ngrid,
                ScalGravModel.SetDensities().origin());
            ScalGravModel.WriteVTK(filename + ".grav.inv.vtk");
            ScalGravModel.WriteNetCDF(filename + ".grav.inv.nc");

          }
        if (HaveScal)
          {
            auto ObsData = ScalGravObjective->GetObservedData();
            jif3D::rvec ScalGravInvData(ScalGravObjective->GetSyntheticData());
            if (ScalGravInvData.size() > 0)
              {
                jif3D::SaveScalarGravityMeasurements(filename + ".inv_sgd.nc",
                    std::vector<double>(ScalGravInvData.begin(), ScalGravInvData.end()),
                    ObsData.GetMeasPosX(), ObsData.GetMeasPosY(), ObsData.GetMeasPosZ(),
                    ScalGravObjective->GetDataError());
                jif3D::Write3DDataToVTK(filename + ".inv_sgd.vtk", "grav_accel",
                    std::vector<double>(ScalGravInvData.begin(), ScalGravInvData.end()),
                    ObsData.GetMeasPosX(), ObsData.GetMeasPosY(), ObsData.GetMeasPosZ());
                jif3D::rvec ScalDiff(ScalGravObjective->GetIndividualMisfit());
                jif3D::SaveScalarGravityMeasurements(filename + ".diff_sgd.nc",
                    std::vector<double>(ScalDiff.begin(), ScalDiff.end()),
                    ObsData.GetMeasPosX(), ObsData.GetMeasPosY(), ObsData.GetMeasPosZ(),
                    ScalGravObjective->GetDataError());
              }
          }
        if (HaveFTG)
          {

            auto ObsData = FTGObjective->GetObservedData();
            jif3D::rvec FTGInvData(FTGObjective->GetSyntheticData());
            if (FTGInvData.size() > 0)
              {
                jif3D::SaveTensorGravityMeasurements(filename + ".inv_ftg.nc",
                    std::vector<double>(FTGInvData.begin(), FTGInvData.end()),
                    ObsData.GetMeasPosX(), ObsData.GetMeasPosY(), ObsData.GetMeasPosZ(),
                    FTGObjective->GetDataError());
                jif3D::Write3DTensorDataToVTK(filename + ".inv_ftg.vtk", "U",
                    std::vector<double>(FTGInvData.begin(), FTGInvData.end()),
                    ObsData.GetMeasPosX(), ObsData.GetMeasPosY(), ObsData.GetMeasPosZ());
              }
          }

      }
  }
