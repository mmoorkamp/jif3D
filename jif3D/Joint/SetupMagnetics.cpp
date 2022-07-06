//============================================================================
// Name        : SetupMagnetics.cpp
// Author      : Mar 1, 2010
// Version     : 
// Copyright   : 2010, mmoorkamp
//============================================================================

#include "SetupMagnetics.h"
#include "../Inversion/ModelTransforms.h"
#include "../Magnetics/ReadWriteMagneticData.h"
#include "../Magnetics/MagneticTransforms.h"
#include "../Magnetics/OMPMagneticSusceptibilityImp.h"
#include "../Gravity/DepthWeighting.h"
#include "../Global/FileUtil.h"
#include "../Global/Noise.h"
#include "../ModelBase/ReadAnyModel.h"
#include "../ModelBase/VTKTools.h"
#include <boost/make_shared.hpp>

namespace jif3D
  {

    const std::string Name = "Susceptibility";
    SetupMagnetics::SetupMagnetics() :
        GeneralDataSetup(Name), inclination(0.0), declination(0.0), fieldstrength(1.0), relerr(
            0.02), minerr(0.0)
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
            "The strength of the inducing magnetic field in nT")("sus_covmod",
            po::value<std::string>(),
            "A file containing the model covariance for susceptibilities")("magdepth",
            "Counteract the decay in sensitivities of magnetic data with depth")("minsus",
            po::value(&minsus)->default_value(-1.0))("maxsus",
            po::value(&maxsus)->default_value(1.0))("logsus",
            "Use a logarithmic transform for susceptibility");

        return desc;
      }

    bool SetupMagnetics::SetupObjective(const po::variables_map &vm,
        jif3D::JointObjective &Objective, jif3D::ThreeDModelBase &InversionMesh,
        jif3D::rvec &CovModVec, std::vector<size_t> &startindices,
        std::vector<std::string> &SegmentNames, std::vector<parametertype> &SegmentTypes,
        boost::filesystem::path TempDir)
      {
        //if we want to use CUDA for forward modeling
        //we set a variable for easier access later and print
        //that information to the screen
        bool wantcuda = false;
        const size_t ngrid = InversionMesh.GetNModelElements();

        if (vm.count("gpu"))
          {
            std::cout << "Using GPU" << "\n";
            wantcuda = true;
          }

        jif3D::TotalFieldMagneticData MagData;

        //we first ask for the weights for scalar and tensor Magnetics
        //as we might not have to do very much if the weights are zero
        std::cout << "Total field magnetics Lambda: ";
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
            std::string magmodelfilename = jif3D::AskFilename(
                "Susceptibility Model Filename: ");
            Model.ReadNetCDF(magmodelfilename);
            if (Model.GetSusceptibilities().num_elements() > 0)
              {
                if (ngrid != Model.GetSusceptibilities().num_elements())
                  {
                    std::string err = "Susceptibility model does not match grid size ! "
                        + std::to_string(Model.GetSusceptibilities().num_elements()) + " "
                        + std::to_string(ngrid);
                    throw jif3D::FatalException(err, __FILE__, __LINE__);
                  }
                StartingParameters.resize(ngrid);
                std::copy(Model.GetSusceptibilities().origin(),
                    Model.GetSusceptibilities().origin() + ngrid,
                    StartingParameters.begin());
              }

            boost::shared_ptr<jif3D::GeneralModelTransform> SusTrans;
            if (vm.count("logsus"))
              {
                jif3D::rvec Ref(ngrid, 1.0);
                SusTrans = boost::make_shared<jif3D::LogTransform>(Ref);
              }
            else
              {
                SusTrans = boost::make_shared<jif3D::TanhTransform>(minsus, maxsus);
              }

            size_t start = startindices.back();
            size_t end = start + ngrid;
            Transform = boost::make_shared<jif3D::MultiSectionTransform>(2 * ngrid, start,
                end, SusTrans);
            startindices.push_back(end);
            SegmentNames.push_back(Name);
            SegmentTypes.push_back(GeneralDataSetup::gridparameter);

            //we want to set the path for temporary file storage
            //the factory function cannot perform this, so we
            //have to assemble the calculator object ourselves
            boost::shared_ptr<
                jif3D::ThreeDGravMagImplementation<jif3D::TotalFieldMagneticData> > Implementation;
            if (wantcuda)
              {
                throw jif3D::FatalException("No GPU support, yet !", __FILE__,
                __LINE__);
              }
            else
              {
                Implementation = boost::shared_ptr<
                    jif3D::ThreeDGravMagImplementation<jif3D::TotalFieldMagneticData> >(
                    new jif3D::OMPMagneticSusceptibilityImp(inclination, declination,
                        fieldstrength));
              }
            CovModVec.resize(ngrid);
            std::fill(CovModVec.begin(), CovModVec.end(), 1.0);
            if (vm.count("sus_covmod"))
              {
                boost::shared_ptr<jif3D::ThreeDModelBase> CovModel = jif3D::ReadAnyModel(
                    vm["sus_covmod"].as<std::string>());

                const size_t ncovmod = CovModel->GetData().num_elements();
                if (ncovmod != ngrid)
                  {
                    throw jif3D::FatalException(
                        " Susceptibility covariance does not have the same number of parameters "
                            + std::to_string(ncovmod) + " as inversion model "
                            + std::to_string(ngrid), __FILE__, __LINE__);

                  }
                else
                  {
                    std::cout << " Setting covariance for susceptibility from file. "
                        << std::endl;
                  }
                std::copy(CovModel->GetData().origin(),
                    CovModel->GetData().origin() + ncovmod, CovModVec.begin());

              }
            if (vm.count("magdepth"))
              {
                boost::shared_ptr<
                    jif3D::ThreeDGravMagImplementation<jif3D::TotalFieldMagneticData> > Implementation =
                    boost::make_shared<jif3D::OMPMagneticSusceptibilityImp>(
                        GetInclination(), GetDeclination(), GetFielStrength());

                jif3D::FullSensitivityGravMagCalculator<jif3D::TotalFieldMagneticData> FullCalc(
                    Implementation);
                FullCalc.SetDataTransform(
                    boost::shared_ptr<jif3D::TotalFieldAnomaly>(
                        new jif3D::TotalFieldAnomaly(GetInclination(),
                            GetDeclination())));
                std::cout << "Calculating depth weighting." << std::endl;
                //now we perform the depth weighting for the sensitivities
                jif3D::rvec SensProfile, WeightVector;
                //we find a measurement site close to the centre of the model and extract the
                //sensitivity variation with depth
                jif3D::rmat Sens;
                jif3D::CalculateMiddleSens(GetModel(), FullCalc, SensProfile);

                double DepthExponent = -3.0;
                //we fit a curve of the form 1/(z+z0)^n to the extracted sensitivities
                double z0 = FitZ0(SensProfile, GetModel().GetZCellSizes(),
                    jif3D::WeightingTerm(DepthExponent));
                std::cout << "Estimated z0: " << z0 << std::endl;
                const size_t zsize = GetModel().GetModelShape()[2];
                //calculate the depth scaling
                jif3D::ConstructDepthWeighting(GetModel().GetZCellSizes(), z0,
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
                jif3D::ThreeDSusceptibilityModel DepthModel(GetModel());
                std::copy(CovModVec.begin(), CovModVec.end(),
                    DepthModel.SetSusceptibilities().origin());
                DepthModel.WriteNetCDF("depth_mag_cov.nc");
                DepthModel.WriteVTK("depth_mag_cov.vtk");
              }
            if (maglambda > JointObjective::MinWeight)
              {
#ifdef MAGDISK
            Calculator = boost::make_shared<MagCalculatorType>(Implementation, TempDir);
#else
                Calculator = boost::make_shared<MagCalculatorType>(Implementation);
#endif
                Calculator->SetDataTransform(
                    boost::shared_ptr<jif3D::TotalFieldAnomaly>(
                        new jif3D::TotalFieldAnomaly(inclination, declination)));

                MagObjective = boost::shared_ptr<
                    jif3D::ThreeDModelObjective<MagCalculatorType> >(
                    new jif3D::ThreeDModelObjective<MagCalculatorType>(*Calculator));
                MagObjective->SetObservedData(MagData);
                MagObjective->SetCoarseModelGeometry(Model);
                std::vector<double> Error(
                    jif3D::ConstructError(MagData.GetData(), MagData.GetErrors(), relerr,
                        minerr));
                MagObjective->SetDataError(Error);

                Objective.AddObjective(MagObjective, Transform, maglambda, "Magnetics",
                    JointObjective::datafit);
                std::cout << " Magnetics ndata: " << MagData.GetData().size()
                    << std::endl;
                std::cout << " Magnetics lambda: " << maglambda << std::endl;
              }
          }

        //indicate whether we added a Magnetics objective function
        //this way the caller can do additional consistency checks
        return (maglambda > JointObjective::MinWeight);
      }

    void SetupMagnetics::IterationOutput(const std::string &filename,
        const jif3D::rvec &ModelVector)
      {
        if (maglambda > JointObjective::MinWeight)
          {
            jif3D::rvec MagInvModel = Transform->GeneralizedToPhysical(ModelVector);
            std::copy(MagInvModel.begin(), MagInvModel.end(),
                Model.SetSusceptibilities().origin());
            Model.WriteVTK(filename + ".mag.inv.vtk");
            Model.WriteNetCDF(filename + ".mag.inv.nc");
          }
      }

    void SetupMagnetics::FinalOutput(const std::string &filename,
        const jif3D::rvec &FinalModelVector)
      {
        if (maglambda > JointObjective::MinWeight)
          {
            jif3D::rvec MagInvModel = Transform->GeneralizedToPhysical(FinalModelVector);
            std::cout << "Writing final susceptibility models " << std::endl;
            std::copy(MagInvModel.begin(), MagInvModel.end(),
                Model.SetSusceptibilities().origin());
            auto ObservedData = GetObjective().GetObservedData();
            jif3D::rvec MagInvData = GetObjective().GetSyntheticData();
            if (MagInvData.size() > 0)
              {
                jif3D::SaveTotalFieldMagneticMeasurements(filename + ".inv_mag.nc",
                    std::vector<double>(MagInvData.begin(), MagInvData.end()),
                    ObservedData.GetMeasPosX(), ObservedData.GetMeasPosY(),
                    ObservedData.GetMeasPosZ(), GetObjective().GetDataError());
                jif3D::Write3DDataToVTK(filename + ".inv_mag.vtk", "T",
                    std::vector<double>(MagInvData.begin(), MagInvData.end()),
                    ObservedData.GetMeasPosX(), ObservedData.GetMeasPosY(),
                    ObservedData.GetMeasPosZ());
                jif3D::rvec MagDiff(GetObjective().GetIndividualMisfit());
                jif3D::SaveTotalFieldMagneticMeasurements(filename + ".diff_mag.nc",
                    std::vector<double>(MagDiff.begin(), MagDiff.end()),
                    ObservedData.GetMeasPosX(), ObservedData.GetMeasPosY(),
                    ObservedData.GetMeasPosZ(), GetObjective().GetDataError());
              }
            Model.WriteVTK(filename + ".mag.inv.vtk");
            Model.WriteNetCDF(filename + ".mag.inv.nc");
          }
      }
  }
