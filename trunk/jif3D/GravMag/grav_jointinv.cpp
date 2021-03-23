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
#include "../ModelBase/ReadAnyModel.h"
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
#include "../Gravity/ScalarGravityData.h"
#include "../Gravity/TensorGravityData.h"
#include "../Magnetics/OMPMagneticImp.h"
#include "../Magnetics/ReadWriteMagneticData.h"
#include "../Magnetics/MagneticTransforms.h"
#include "../Magnetics/MagneticData.h"
#include "../MI/MutualInformationConstraint.h"
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

    double mindens, maxdens, minsus, maxsus;
    int mibins = 10;
    po::options_description desc("General options");
    desc.add_options()("help", "produce help message")("threads", po::value<int>(),
        "The number of openmp threads")("dens_covmod", po::value<std::string>(),
        "A file containing the model covariance")("magdepth",
        "Counteract the decay in sensitivities of magnetic data with depth")("gravdepth",
        "Counteract the decay in sensitivities of gravity data with depth")("mindens",
        po::value(&mindens)->default_value(-1000.0))("maxdens",
        po::value(&maxdens)->default_value(1000.0))("minsus",
        po::value(&minsus)->default_value(-1.0))("maxsus",
        po::value(&maxsus)->default_value(1.0))("mutual_information", po::value(&mibins),
        "Use mutual information coupling, specify number of bins in histogram");

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
    boost::shared_ptr<jif3D::ThreeDModelBase> Mesh = jif3D::ReadAnyModel(meshfilename);
    const size_t ngrid = Mesh->GetNModelElements();

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

    boost::shared_ptr<jif3D::GeneralModelTransform> Copy = boost::make_shared<
        jif3D::ModelCopyTransform>();

    jif3D::rvec InvModel(2 * ngrid);

    boost::shared_ptr<jif3D::GeneralModelTransform> GravityRegTransform =
        boost::make_shared<jif3D::MultiSectionTransform>(2 * ngrid, 0, ngrid, Copy);
    boost::shared_ptr<jif3D::GeneralModelTransform> MagneticsRegTransform =
        boost::make_shared<jif3D::MultiSectionTransform>(2 * ngrid, ngrid, 2 * ngrid,
            Copy);

    boost::shared_ptr<jif3D::GeneralModelTransform> DensTrans = boost::make_shared<
        jif3D::TanhTransform>(mindens, maxdens);
    boost::shared_ptr<jif3D::GeneralModelTransform> SusTrans = boost::make_shared<
        jif3D::TanhTransform>(minsus, maxsus);
    boost::shared_ptr<jif3D::GeneralModelTransform> GravityTransform = boost::make_shared<
        jif3D::MultiSectionTransform>(2 * ngrid, 0, ngrid, DensTrans);
    boost::shared_ptr<jif3D::GeneralModelTransform> MagneticsTransform =
        boost::make_shared<jif3D::MultiSectionTransform>(2 * ngrid, ngrid, 2 * ngrid,
            SusTrans);

    auto Objective = boost::make_shared<jif3D::JointObjective>(true);
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

    boost::shared_ptr<jif3D::CrossGradient> GravMagCross(new jif3D::CrossGradient(*Mesh));
    boost::shared_ptr<jif3D::MultiSectionTransform> GravMagTrans(
        new jif3D::MultiSectionTransform(2 * ngrid));
    GravMagTrans->AddSection(0, ngrid, Copy);
    GravMagTrans->AddSection(ngrid, 2 * ngrid, Copy);

    if (HaveGrav)
      {
        jif3D::rvec GravMod(GravitySetup.GetScalModel().GetDensities().num_elements());
        std::copy(GravitySetup.GetScalModel().GetDensities().origin(),
            GravitySetup.GetScalModel().GetDensities().origin() + ngrid, GravMod.begin());
        jif3D::rvec GravGen = DensTrans->PhysicalToGeneralized(GravMod);
        std::copy(GravGen.begin(), GravGen.end(), InvModel.begin());
      }
    else
      {
        std::fill_n(InvModel.begin(), ngrid, 0.0);
      }

    if (HaveMag)
      {
        jif3D::rvec MagMod(
            MagneticsSetup.GetModel().GetSusceptibilities().num_elements());
        std::copy(MagneticsSetup.GetModel().GetSusceptibilities().origin(),
            MagneticsSetup.GetModel().GetSusceptibilities().origin() + ngrid,
            MagMod.begin());
        jif3D::rvec MagGen = SusTrans->PhysicalToGeneralized(MagMod);
        std::copy(MagGen.begin(), MagGen.end(), InvModel.begin() + ngrid);

      }
    else
      {
        std::fill_n(InvModel.begin() + ngrid, ngrid, 0.0);
      }

    boost::shared_ptr<jif3D::ObjectiveFunction> Regularization = RegSetup.SetupObjective(
        vm, *Mesh, CovModVec);

    CovModVec.resize(2 * ngrid, true);
    std::fill(CovModVec.begin(), CovModVec.end(), 1.0);

    if (vm.count("magdepth"))
      {
        boost::shared_ptr<jif3D::ThreeDGravMagImplementation<jif3D::MagneticData> > Implementation(
            new jif3D::OMPMagneticImp(MagneticsSetup.GetInclination(),
                MagneticsSetup.GetDeclination(), MagneticsSetup.GetFielStrength()));
        jif3D::FullSensitivityGravMagCalculator<jif3D::MagneticData> FullCalc(
            Implementation);
        FullCalc.SetDataTransform(
            boost::shared_ptr<jif3D::TotalFieldAnomaly>(
                new jif3D::TotalFieldAnomaly(MagneticsSetup.GetInclination(),
                    MagneticsSetup.GetDeclination())));
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

        boost::shared_ptr<jif3D::ThreeDGravMagImplementation<jif3D::ScalarGravityData> > Implementation(
            new jif3D::ScalarOMPGravityImp);
        jif3D::FullSensitivityGravMagCalculator<jif3D::ScalarGravityData> FullCalc(
            Implementation);

        std::cout << "Calculating depth weighting." << std::endl;
        //now we perform the depth weighting for the sensitivities
        jif3D::rvec SensProfile, WeightVector;
        //we find a measurement site close to the centre of the model and extract the
        //sensitivity variation with depth
        jif3D::rmat Sens;
        jif3D::CalculateMiddleSens(GravitySetup.GetScalModel(), FullCalc, SensProfile);

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
        Objective->AddObjective(Regularization, GravityRegTransform, gravreglambda,
            "GravReg", jif3D::JointObjective::regularization);
      }
    double magreglambda = 1.0;
    std::cout << "Magnetics Regularization Lambda: ";
    std::cin >> magreglambda;
    if (magreglambda > 0.0)
      {
        Objective->AddObjective(Regularization, MagneticsRegTransform, magreglambda,
            "MagReg", jif3D::JointObjective::regularization);
      }

    if (vm.count("mutual_information"))
      {
        boost::shared_ptr<jif3D::MutualInformationConstraint> GravMagMI =
            boost::make_shared<jif3D::MutualInformationConstraint>(-2.0, 2.0, -2.0, 2.0,
                mibins);
        double milambda = 10.0;
        std::cout << "MI weight: ";
        std::cin >> milambda;
        if (milambda > 0.0)
          {
            Objective->AddObjective(GravMagMI, GravMagTrans, milambda, "MI",
                jif3D::JointObjective::coupling);
          }
      }
    else
      {
        double crosslambda = 10.0;
        std::cout << "Cross-gradient weight: ";
        std::cin >> crosslambda;
        if (crosslambda > 0.0)
          {
            Objective->AddObjective(GravMagCross, GravMagTrans, crosslambda, "Cross",
                jif3D::JointObjective::coupling);
          }
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
    std::ofstream rmsfile("rms.out");

    double InitialMisfit = Objective->CalcMisfit(InvModel);
    StoreMisfit(misfitfile, 0, InitialMisfit, *Objective);

    jif3D::ThreeDGravityModel GravModel(GravitySetup.GetScalModel());
    jif3D::ThreeDMagneticModel MagModel(MagneticsSetup.GetModel());
    bool terminate = false;
    while (iteration < maxiter && !terminate)
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

            StoreMisfit(misfitfile, iteration, Optimizer->GetMisfit(), *Objective);
            StoreRMS(rmsfile, iteration, *Objective);
            std::cout << "\n\n";
          } catch (jif3D::FatalException &e)
          {
            std::cerr << e.what() << std::endl;
            iteration = maxiter;
          }
        //we stop when either we do not make any improvement any more
        terminate = CheckConvergence(*Objective);
        //or the file abort exists in the current directory
        terminate = terminate || jif3D::WantAbort();
      }

    jif3D::rvec DensInvModel(GravityTransform->GeneralizedToPhysical(InvModel));
    jif3D::rvec MagInvModel(MagneticsTransform->GeneralizedToPhysical(InvModel));

    //calculate the predicted data
    std::cout << "Calculating response of inversion model." << std::endl;
    if (GravitySetup.GetHaveScal())
      {
        auto ObsData = GravitySetup.GetScalGravObjective().GetObservedData();
        jif3D::rvec ScalGravInvData(
            GravitySetup.GetScalGravObjective().GetSyntheticData());

        jif3D::SaveScalarGravityMeasurements(modelfilename + ".inv_sgd.nc",
            std::vector<double>(ScalGravInvData.begin(), ScalGravInvData.end()),
            ObsData.GetMeasPosX(), ObsData.GetMeasPosY(), ObsData.GetMeasPosZ(),
            GravitySetup.GetScalGravObjective().GetDataError());
        jif3D::Write3DDataToVTK(modelfilename + ".inv_sgd.vtk", "grav_accel",
            std::vector<double>(ScalGravInvData.begin(), ScalGravInvData.end()),
            ObsData.GetMeasPosX(), ObsData.GetMeasPosY(), ObsData.GetMeasPosZ());
        jif3D::rvec ScalDiff(GravitySetup.GetScalGravObjective().GetIndividualMisfit());
        jif3D::SaveScalarGravityMeasurements(modelfilename + ".diff_sgd.nc",
            std::vector<double>(ScalDiff.begin(), ScalDiff.end()), ObsData.GetMeasPosX(),
            ObsData.GetMeasPosY(), ObsData.GetMeasPosZ(),
            GravitySetup.GetScalGravObjective().GetDataError());
      }
    if (GravitySetup.GetHaveFTG())
      {
        std::copy(DensInvModel.begin(), DensInvModel.begin() + ngrid,
            GravModel.SetDensities().origin());
        auto ObsData = GravitySetup.GetFTGObjective().GetObservedData();
        jif3D::rvec FTGInvData(GravitySetup.GetFTGObjective().GetSyntheticData());
        jif3D::SaveTensorGravityMeasurements(modelfilename + ".inv_ftg.nc",
            std::vector<double>(FTGInvData.begin(), FTGInvData.end()),
            ObsData.GetMeasPosX(), ObsData.GetMeasPosY(), ObsData.GetMeasPosZ(),
            GravitySetup.GetScalGravObjective().GetDataError());
        jif3D::Write3DTensorDataToVTK(modelfilename + ".inv_ftg.vtk", "U",
            std::vector<double>(FTGInvData.begin(), FTGInvData.end()),
            ObsData.GetMeasPosX(), ObsData.GetMeasPosY(), ObsData.GetMeasPosZ());
      }

    if (HaveMag)
      {
        std::copy(MagInvModel.begin(), MagInvModel.begin() + ngrid,
            MagModel.SetSusceptibilities().origin());
        auto ObservedData = MagneticsSetup.GetObjective().GetObservedData();
        jif3D::rvec MagInvData = MagneticsSetup.GetObjective().GetSyntheticData();
        jif3D::SaveTotalFieldMagneticMeasurements(modelfilename + ".inv_mag.nc",
            std::vector<double>(MagInvData.begin(), MagInvData.end()),
            ObservedData.GetMeasPosX(), ObservedData.GetMeasPosY(),
            ObservedData.GetMeasPosZ(), MagneticsSetup.GetObjective().GetDataError());
        jif3D::Write3DDataToVTK(modelfilename + ".inv_mag.vtk", "T",
            std::vector<double>(MagInvData.begin(), MagInvData.end()),
            ObservedData.GetMeasPosX(), ObservedData.GetMeasPosY(),
            ObservedData.GetMeasPosZ());
        jif3D::rvec MagDiff(MagneticsSetup.GetObjective().GetIndividualMisfit());
        jif3D::SaveTotalFieldMagneticMeasurements(modelfilename + ".diff_mag.nc",
            std::vector<double>(MagDiff.begin(), MagDiff.end()),
            ObservedData.GetMeasPosX(), ObservedData.GetMeasPosY(),
            ObservedData.GetMeasPosZ(), MagneticsSetup.GetObjective().GetDataError());
        MagModel.WriteVTK(modelfilename + ".mag.inv.vtk");
        MagModel.WriteNetCDF(modelfilename + ".mag.inv.nc");
      }

    if (vm.count("mutual_information"))
      {

      }
    else
      {

        auto ObjectiveTypes = Objective->GetObjectiveTypes();
        std::vector<std::string> Names = Objective->GetObjectiveNames();
        for (size_t i = 0; i < ObjectiveTypes.size(); ++i)
          {
            if (ObjectiveTypes.at(i) == jif3D::JointObjective::coupling)
              {
                jif3D::rvec CG(Objective->GetObjective(i).GetIndividualMisfit());
                const size_t nx = Mesh->GetData().shape()[0];
                const size_t ny = Mesh->GetData().shape()[1];
                const size_t nz = Mesh->GetData().shape()[2];
                const size_t nmod = nx * ny * nz;
                jif3D::ThreeDModelBase::t3DModelData XGrad(boost::extents[nx][ny][nz]);
                jif3D::ThreeDModelBase::t3DModelData YGrad(boost::extents[nx][ny][nz]);
                jif3D::ThreeDModelBase::t3DModelData ZGrad(boost::extents[nx][ny][nz]);
                std::copy(CG.begin(), CG.begin() + nmod, XGrad.origin());
                std::copy(CG.begin() + nmod, CG.begin() + 2 * nmod, YGrad.origin());
                std::copy(CG.begin() + 2 * nmod, CG.begin() + 3 * nmod, ZGrad.origin());
                std::string Name = Names.at(i);
                jif3D::Write3DVectorModelToVTK(Name + ".vtk", Name,
                    Mesh->GetXCoordinates(), Mesh->GetYCoordinates(),
                    Mesh->GetZCoordinates(), XGrad, YGrad, ZGrad);
                jif3D::ThreeDGravityModel CGModel(GravitySetup.GetScalModel());
                jif3D::ThreeDModelBase::t3DModelData AbsCG(boost::extents[nx][ny][nz]);

                std::transform(CG.begin(), CG.begin() + nmod, AbsCG.origin(),
                    [](double val)
                      { return val*val;});
                std::transform(CG.begin() + nmod, CG.begin() + 2 * nmod, AbsCG.origin(),
                    AbsCG.origin(), [](double val1, double val2)
                      { return val2 + val1*val1;});
                std::transform(CG.begin() + 2 * nmod, CG.begin() + 3 * nmod,
                    AbsCG.origin(), AbsCG.origin(), [](double val1, double val2)
                      { return val2 + val1*val1;});
                std::transform(AbsCG.origin(), AbsCG.origin() + nmod, AbsCG.origin(),
                    [](double val)
                      { return std::sqrt(val);});
                CGModel.SetDensities() = AbsCG;
                CGModel.WriteNetCDF(Name+".cg_abs.nc");
                CGModel.WriteVTK(Name+".cg_abs.vtk");
                std::vector<double> CG_Cov(Objective->GetObjective(i).GetDataError());
                std::copy(CG_Cov.begin(), CG_Cov.begin() + nmod, XGrad.origin());
                std::copy(CG_Cov.begin() + nmod, CG_Cov.begin() + 2 * nmod,
                    YGrad.origin());
                std::copy(CG_Cov.begin() + 2 * nmod, CG_Cov.begin() + 3 * nmod,
                    ZGrad.origin());
                jif3D::Write3DVectorModelToVTK(Name + "_cov.vtk", Name,
                    Mesh->GetXCoordinates(), Mesh->GetYCoordinates(),
                    Mesh->GetZCoordinates(), XGrad, YGrad, ZGrad);
              }
          }

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
    std::cout << std::endl;
  }
