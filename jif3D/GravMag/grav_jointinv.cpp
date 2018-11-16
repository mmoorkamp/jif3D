//============================================================================
// Name        : gravinv.cpp
// Author      : Max Moorkamp
// Version     :
// Copyright   : 2008, MM
//============================================================================

#include "../Global/Serialization.h"
#include "../Global/convert.h"
#include "../Global/FatalException.h"
#include "../Global/NumUtil.h"
#include "../Global/VectorTransform.h"
#include "../Global/FileUtil.h"
#include "../Global/Noise.h"
#include "../ModelBase/VTKTools.h"
#include "../ModelBase/NetCDFModelTools.h"
#include "../Inversion/LimitedMemoryQuasiNewton.h"
#include "../Inversion/NonLinearConjugateGradient.h"
#include "../Inversion/JointObjective.h"
#include "../Inversion/DiagonalCovariance.h"
#include "../Regularization/MinDiffRegularization.h"
#include "../Regularization/CrossGradient.h"
#include "../Inversion/ModelTransforms.h"
#include "../Gravity/ReadWriteGravityData.h"
#include "../GravMag/ThreeDGravMagCalculator.h"
#include "../GravMag/MinMemGravMagCalculator.h"
#include "../Gravity/DepthWeighting.h"
#include "../Gravity/ThreeDGravityFactory.h"
#include "../Magnetics/OMPMagneticImp.h"
#include "../Magnetics/ReadWriteMagneticData.h"
#include "../Magnetics/MagneticTransforms.h"
#include "../Joint/SetupRegularization.h"
#include "../Joint/SetupInversion.h"
#include "../Joint/SetupGravity.h"
#include "../Joint/SetupMagnetics.h"
#include "../Joint/InversionOutput.h"
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#ifdef HAVEOPENMP
#include <omp.h>
#endif
#include <boost/program_options.hpp>
namespace ublas = boost::numeric::ublas;
namespace po = boost::program_options;

int main(int argc, char *argv[])
  {

    po::options_description desc("General options");
    desc.add_options()("help", "produce help message")("threads", po::value<int>(),
        "The number of openmp threads")("dens_covmod", po::value<std::string>(),
        "A file containing the model covariance")("magdepth",
        "Counteract the decay in sensitivities of magnetic data with depth")("gravdepth",
        "Counteract the decay in sensitivities of gravity data with depth");

    jif3D::SetupRegularization RegSetup;
    jif3D::SetupInversion InversionSetup;
    jif3D::SetupGravity GravitySetup;
    jif3D::SetupMagnetics MagneticsSetup;
    desc.add(RegSetup.SetupOptions());
    desc.add(InversionSetup.SetupOptions());
    desc.add(GravitySetup.SetupOptions());
    desc.add(MagneticsSetup.SetupOptions());
    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);

    if (vm.count("help"))
      {
        std::cout << desc << "\n";
        return 1;
      }
#ifdef HAVEOPENMP
    if (vm.count("threads"))
      {
        omp_set_num_threads(vm["threads"].as<int>());
      }
#endif

    std::string meshfilename = jif3D::AskFilename("Mesh filename: ");
    jif3D::ThreeDGravityModel Mesh;
    Mesh.ReadNetCDF(meshfilename);
    const size_t ngrid = Mesh.GetNModelElements();

    jif3D::rvec CovModVec;
    if (vm.count("dens_covmod"))
      {
        jif3D::ThreeDGravityModel CovModel;
        CovModel.ReadNetCDF(vm["dens_covmod"].as<std::string>());
        const size_t ncovmod = CovModel.GetDensities().num_elements();
        std::copy(CovModel.GetDensities().origin(),
            CovModel.GetDensities().origin() + ncovmod, CovModVec.begin());
      }

    jif3D::rvec RefVec(ngrid);
    std::fill(RefVec.begin(), RefVec.end(), 1.0);

    boost::shared_ptr<jif3D::GeneralModelTransform> Trans = boost::shared_ptr<
        jif3D::GeneralModelTransform>(new jif3D::NormalizeTransform(RefVec));
    jif3D::rvec InvModel(2 * ngrid);

    boost::shared_ptr<jif3D::GeneralModelTransform> GravityTransform(
        new jif3D::MultiSectionTransform(2 * ngrid, 0, ngrid, Trans));
    boost::shared_ptr<jif3D::GeneralModelTransform> MagneticsTransform(
        new jif3D::MultiSectionTransform(2 * ngrid, ngrid, 2 * ngrid, Trans));

    boost::shared_ptr<jif3D::JointObjective> Objective(new jif3D::JointObjective(true));
    bool HaveGrav = GravitySetup.SetupObjective(vm, *Objective.get(), GravityTransform);
    bool HaveMag = MagneticsSetup.SetupObjective(vm, *Objective.get(),
        MagneticsTransform);

    if (HaveMag)
      {
        if (ngrid != MagneticsSetup.GetModel().GetSusceptibilities().num_elements())
          {
            std::cerr << "Gravity model and magnetic model have different size ! "
                << std::endl;
            return 100;
          }
      }

    boost::shared_ptr<jif3D::CrossGradient> GravMagCross(new jif3D::CrossGradient(Mesh));
    boost::shared_ptr<jif3D::MultiSectionTransform> GravMagTrans(
        new jif3D::MultiSectionTransform(2 * ngrid));
    GravMagTrans->AddSection(0, ngrid, Trans);
    GravMagTrans->AddSection(ngrid, 2 * ngrid, Trans);

    if (HaveGrav)
      {
        std::copy(GravitySetup.GetScalModel().GetDensities().origin(),
            GravitySetup.GetScalModel().GetDensities().origin() + ngrid,
            InvModel.begin());
      }
    else
      {
        std::fill_n(InvModel.begin(), ngrid, 0.0);
      }

    if (HaveMag)
      {
        std::copy(MagneticsSetup.GetModel().GetSusceptibilities().origin(),
            MagneticsSetup.GetModel().GetSusceptibilities().origin() + ngrid,
            InvModel.begin() + ngrid);
      }
    else
      {
        std::fill_n(InvModel.begin() + ngrid, ngrid, 0.0);
      }

    boost::shared_ptr<jif3D::ObjectiveFunction> Regularization = RegSetup.SetupObjective(
        vm, Mesh, CovModVec);

    CovModVec.resize(2 * ngrid, true);
    std::fill(CovModVec.begin(), CovModVec.end(), 1.0);

    if (vm.count("magdepth"))
      {
        boost::shared_ptr<jif3D::ThreeDGravMagImplementation<jif3D::ThreeDMagneticModel> > Implementation(
            new jif3D::OMPMagneticImp(MagneticsSetup.GetInclination(),
                MagneticsSetup.GetDeclination(), MagneticsSetup.GetFielStrength()));
        jif3D::FullSensitivityGravMagCalculator<jif3D::ThreeDMagneticModel> FullCalc(
            Implementation);
        FullCalc.SetDataTransform(
            boost::shared_ptr<jif3D::TotalFieldAnomaly>(
                new jif3D::TotalFieldAnomaly(MagneticsSetup.GetInclination(),
                    MagneticsSetup.GetDeclination(), MagneticsSetup.GetFielStrength())));
        std::cout << "Calculating depth weighting." << std::endl;
        //now we perform the depth weighting for the sensitivities
        jif3D::rvec SensProfile, WeightVector;
        //we find a measurement site close to the centre of the model and extract the
        //sensitivity variation with depth
        jif3D::rmat Sens;
        jif3D::CalculateMiddleSens(MagneticsSetup.GetModel(), FullCalc, SensProfile);

        double DepthExponent = -3.0;
        //we fit a curve of the form 1/(z+z0)^n to the extracted sensitivities
        double z0 = FitZ0(SensProfile, MagneticsSetup.GetModel().GetZCellSizes(),
            jif3D::WeightingTerm(DepthExponent));
        std::cout << "Estimated z0: " << z0 << std::endl;
        const size_t zsize = MagneticsSetup.GetModel().GetModelShape()[2];
        //calculate the depth scaling
        jif3D::ConstructDepthWeighting(MagneticsSetup.GetModel().GetZCellSizes(), z0,
            WeightVector, jif3D::WeightingTerm(DepthExponent));
        for (size_t i = 0; i < ngrid; ++i)
          {
            CovModVec(ngrid + i) = WeightVector(i % zsize);
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
      }

    if (vm.count("gravdepth"))
      {

        boost::shared_ptr<jif3D::ThreeDGravMagImplementation<jif3D::ThreeDGravityModel> > Implementation(
            new jif3D::ScalarOMPGravityImp);
        jif3D::FullSensitivityGravMagCalculator<jif3D::ThreeDGravityModel> FullCalc(
            Implementation);

        std::cout << "Calculating depth weighting." << std::endl;
        //now we perform the depth weighting for the sensitivities
        jif3D::rvec SensProfile, WeightVector;
        //we find a measurement site close to the centre of the model and extract the
        //sensitivity variation with depth
        jif3D::rmat Sens;
        jif3D::CalculateMiddleSens(GravitySetup.GetScalModel(), FullCalc,
            SensProfile);

        double DepthExponent = -2.0;
        //we fit a curve of the form 1/(z+z0)^n to the extracted sensitivities
        double z0 = FitZ0(SensProfile, GravitySetup.GetScalModel().GetZCellSizes(),
            jif3D::WeightingTerm(DepthExponent));
        std::cout << "Estimated z0: " << z0 << std::endl;
        const size_t zsize = GravitySetup.GetScalModel().GetModelShape()[2];
        //calculate the depth scaling
        jif3D::ConstructDepthWeighting(GravitySetup.GetScalModel().GetZCellSizes(), z0,
            WeightVector, jif3D::WeightingTerm(DepthExponent));
        for (size_t i = 0; i < ngrid; ++i)
          {
            CovModVec(i) = WeightVector(i % zsize);
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
      }

    double gravreglambda = 1.0;
    std::cout << "Gravity Regularization Lambda: ";
    std::cin >> gravreglambda;
    if (gravreglambda > 0.0)
      {
        Objective->AddObjective(Regularization, GravityTransform, gravreglambda,
            "GravReg", jif3D::JointObjective::regularization);
      }
    double magreglambda = 1.0;
    std::cout << "Magnetics Regularization Lambda: ";
    std::cin >> magreglambda;
    if (magreglambda > 0.0)
      {
        Objective->AddObjective(Regularization, MagneticsTransform, magreglambda,
            "MagReg", jif3D::JointObjective::regularization);
      }

    double crosslambda = 10.0;
    std::cout << "Cross-gradient weight: ";
    std::cin >> crosslambda;
    if (crosslambda > 0.0)
      {
        Objective->AddObjective(GravMagCross, GravMagTrans, crosslambda, "Cross",
            jif3D::JointObjective::coupling);
      }
    std::cout << "Performing inversion." << std::endl;

    //jif3D::rvec Ones(CovModVec);
    //std::fill(Ones.begin(),Ones.end(),1.0);
    auto CovObj = boost::make_shared<jif3D::DiagonalCovariance>(CovModVec);

    boost::shared_ptr<jif3D::GradientBasedOptimization> Optimizer =
        InversionSetup.ConfigureInversion(vm, Objective, InvModel, CovObj);

    size_t iteration = 0;
    size_t maxiter = 30;
    std::cout << "Maximum number of iterations: ";
    std::cin >> maxiter;
    std::string modelfilename("result");
    std::ofstream misfitfile("misfit.out");

    double InitialMisfit = Objective->CalcMisfit(InvModel);
    StoreMisfit(misfitfile, 0, InitialMisfit, *Objective);

    jif3D::ThreeDGravityModel GravModel(GravitySetup.GetScalModel());
    jif3D::ThreeDMagneticModel MagModel(MagneticsSetup.GetModel());
    bool terminate = false;
    do
      {
        try
          {
            std::cout << "Iteration: " << iteration << std::endl;
            Optimizer->MakeStep(InvModel);

            if (HaveGrav)
              {
                jif3D::rvec DensInvModel = GravityTransform->GeneralizedToPhysical(
                    InvModel);
                std::copy(DensInvModel.begin(), DensInvModel.end(),
                    GravModel.SetDensities().origin());
                GravModel.WriteVTK(
                    modelfilename + jif3D::stringify(iteration) + ".grav.inv.vtk");
                GravModel.WriteNetCDF(
                    modelfilename + jif3D::stringify(iteration) + ".grav.inv.nc");
              }
            if (HaveMag)
              {
                jif3D::rvec MagInvModel = MagneticsTransform->GeneralizedToPhysical(
                    InvModel);
                std::copy(MagInvModel.begin(), MagInvModel.end(),
                    MagModel.SetSusceptibilities().origin());
                MagModel.WriteVTK(
                    modelfilename + jif3D::stringify(iteration) + ".mag.inv.vtk");
                MagModel.WriteNetCDF(
                    modelfilename + jif3D::stringify(iteration) + ".mag.inv.nc");
              }
            ++iteration;
            std::cout << "Gradient Norm: " << Optimizer->GetGradNorm() << std::endl;
            std::cout << "Currrent Misfit: " << Optimizer->GetMisfit() << std::endl;
            std::cout << "Currrent Gradient: " << Optimizer->GetGradNorm() << std::endl;
            misfitfile << std::setw(5) << iteration << " " << std::setw(15)
                << Optimizer->GetMisfit() << " ";
            for (size_t i = 0; i < Objective->GetIndividualFits().size(); ++i)
              {
                misfitfile << std::setw(15) << Objective->GetIndividualFits().at(i)
                    << " ";
              }

            misfitfile << " " << Objective->GetNEval();
            misfitfile << std::endl;
          } catch (jif3D::FatalException &e)
          {
            std::cerr << e.what() << std::endl;
            iteration = maxiter;
          }
        //we stop when either we do not make any improvement any more
        terminate = CheckConvergence(*Objective);
        //or the file abort exists in the current directory
        terminate = terminate || jif3D::WantAbort();
      } while (iteration < maxiter && !terminate);

    jif3D::rvec DensInvModel(GravityTransform->GeneralizedToPhysical(InvModel));
    jif3D::rvec MagInvModel(MagneticsTransform->GeneralizedToPhysical(InvModel));

    //calculate the predicted data
    std::cout << "Calculating response of inversion model." << std::endl;
    typedef typename jif3D::MinMemGravMagCalculator<jif3D::ThreeDGravityModel> GravCalculatorType;
    typedef typename jif3D::MinMemGravMagCalculator<jif3D::ThreeDMagneticModel> MagCalculatorType;
    if (GravitySetup.GetHaveScal())
      {
        std::copy(DensInvModel.begin(), DensInvModel.begin() + ngrid,
            GravModel.SetDensities().origin());
        boost::shared_ptr<GravCalculatorType> ScalGravityCalculator = boost::shared_ptr<
            GravCalculatorType>(
            jif3D::CreateGravityCalculator<GravCalculatorType>::MakeScalar());
        jif3D::rvec GravInvData(ScalGravityCalculator->Calculate(GravModel));
        jif3D::SaveScalarGravityMeasurements(modelfilename + ".inv_sgd.nc", GravInvData,
            GravModel.GetMeasPosX(), GravModel.GetMeasPosY(), GravModel.GetMeasPosZ(),
            GravitySetup.GetScalGravObjective().GetDataError());
        jif3D::Write3DDataToVTK(modelfilename + ".inv_sgd.vtk", "grav_accel", GravInvData,
            GravModel.GetMeasPosX(), GravModel.GetMeasPosY(), GravModel.GetMeasPosZ());
      }
    if (GravitySetup.GetHaveFTG())
      {
        std::copy(DensInvModel.begin(), DensInvModel.begin() + ngrid,
            GravModel.SetDensities().origin());
        boost::shared_ptr<GravCalculatorType> FTGGravityCalculator = boost::shared_ptr<
            GravCalculatorType>(
            jif3D::CreateGravityCalculator<GravCalculatorType>::MakeTensor());
        jif3D::rvec FTGInvData(FTGGravityCalculator->Calculate(GravModel));
        jif3D::SaveTensorGravityMeasurements(modelfilename + ".inv_ftg.nc", FTGInvData,
            GravModel.GetMeasPosX(), GravModel.GetMeasPosY(), GravModel.GetMeasPosZ(),
            GravitySetup.GetFTGObjective().GetDataError());
        jif3D::Write3DTensorDataToVTK(modelfilename + ".inv_ftg.vtk", "U", FTGInvData,
            GravModel.GetMeasPosX(), GravModel.GetMeasPosY(), GravModel.GetMeasPosZ());
      }

    if (HaveMag)
      {
        std::copy(MagInvModel.begin(), MagInvModel.begin() + ngrid,
            MagModel.SetSusceptibilities().origin());

        jif3D::rvec MagInvData(MagneticsSetup.GetCalculator()->Calculate(MagModel));
        jif3D::SaveTotalFieldMagneticMeasurements(modelfilename + ".inv_mag.nc",
            MagInvData, MagModel.GetMeasPosX(), MagModel.GetMeasPosY(),
            MagModel.GetMeasPosZ(), MagneticsSetup.GetObjective().GetDataError());
        jif3D::Write3DDataToVTK(modelfilename + ".inv_mag.vtk", "T", MagInvData,
            MagModel.GetMeasPosX(), MagModel.GetMeasPosY(), MagModel.GetMeasPosZ());
        MagModel.WriteVTK(modelfilename + ".mag.inv.vtk");
        MagModel.WriteNetCDF(modelfilename + ".mag.inv.nc");
      }
    //and write out the data and model
    //here we have to distinguish again between scalar and ftg data
    std::cout << "Writing out inversion results." << std::endl;
    if (HaveGrav)
      {
        GravModel.WriteVTK(modelfilename + ".grav.inv.vtk");
        GravModel.WriteNetCDF(modelfilename + ".grav.inv.nc");
      }
    if (HaveMag)
      {
        MagModel.WriteVTK(modelfilename + ".mag.inv.vtk");
        MagModel.WriteNetCDF(modelfilename + ".mag.inv.nc");
      }
    std::cout << std::endl;
  }
