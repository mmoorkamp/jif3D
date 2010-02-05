//============================================================================
// Name        : jointinv.cpp
// Author      : May 12, 2009
// Version     :
// Copyright   : 2009, mmoorkamp
//============================================================================

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <omp.h>
#include <boost/bind.hpp>
#include <boost/program_options.hpp>
#include "../Global/convert.h"
#include "../Global/FatalException.h"
#include "../Global/NumUtil.h"
#include "../Global/VectorTransform.h"
#include "../Global/FileUtil.h"
#include "../Global/Noise.h"
#include "../ModelBase/VTKTools.h"
#include "../ModelBase/NetCDFTools.h"
#include "../Inversion/LimitedMemoryQuasiNewton.h"
#include "../Inversion/NonLinearConjugateGradient.h"
#include "../Inversion/JointObjective.h"
#include "../Regularization/MinDiffRegularization.h"
#include "../Regularization/GradientRegularization.h"
#include "../Regularization/CurvatureRegularization.h"
#include "../Inversion/ModelTransforms.h"
#include "../Tomo/ThreeDSeismicModel.h"
#include "../Tomo/ReadWriteTomographyData.h"
#include "../Tomo/TomographyObjective.h"
#include "../Tomo/TomographyCalculator.h"
#include "../Gravity/GravityObjective.h"
#include "../Gravity/ReadWriteGravityData.h"
#include "../Gravity/ThreeDGravityCalculator.h"
#include "../Gravity/MinMemGravityCalculator.h"
#include "../Gravity/DepthWeighting.h"
#include "../Gravity/ThreeDGravityFactory.h"
#include "../MT/X3DObjective.h"
#include "../MT/X3DModel.h"
#include "../MT/ReadWriteImpedances.h"
#include <boost/date_time/posix_time/posix_time.hpp>

namespace ublas = boost::numeric::ublas;
namespace po = boost::program_options;

int main(int argc, char *argv[])
  {

    po::options_description desc("Allowed options");
    desc.add_options()("help", "produce help message")("cpu",
        "Perform calculation on CPU [default]")("gpu",
        "Perform calculation on GPU")("threads", po::value<int>(),
        "The number of openmp threads")("corrpairs", po::value<int>(),
        "The number correction pairs for L-BFGS")("nlcg",
        "Use NLCG optimization")("xreg", po::value<double>(),
        "The weight for the regularization in x-direction")("yreg", po::value<
        double>(), "The weight for the regularization in y-direction")("zreg",
        po::value<double>(), "The weight for the regularization in z-direction")(
        "curvreg", "Use model curvature for regularization.")("refmod",
        "Substract the starting model as a reference model for the regularization.");

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);

    if (vm.count("help"))
      {
        std::cout << desc << "\n";
        return 1;
      }
    bool wantcuda = false;

    if (vm.count("gpu"))
      {
        std::cout << "Using GPU" << "\n";
        wantcuda = true;
      }
    if (vm.count("threads"))
      {
        omp_set_num_threads(vm["threads"].as<int> ());
      }
    int correctionpairs = 5;
    if (vm.count("corrpairs"))
      {
        correctionpairs = vm["corrpairs"].as<int> ();
      }
    bool wantnlcg = false;
    if (vm.count("nlcg"))
      {
        wantnlcg = true;
      }
    bool wantcurvreg = false;
    if (vm.count("curvreg"))
      {
        wantcurvreg = true;
      }

    //these objects hold information about the measurements and their geometry
    jiba::rvec TomoData, ScalGravData, FTGData, MTData;

    //first we read in the starting model and the measured data
    std::string modelfilename = jiba::AskFilename("Starting model Filename: ");
    //we read in the starting modelfile
    jiba::ThreeDSeismicModel TomoModel, TomoFineGeometry;
    TomoModel.ReadNetCDF(modelfilename);
    TomoModel.WriteVTK(modelfilename + ".vtk");
    std::string tomogeometryfilename = jiba::AskFilename("Tomography forward geometry filename: ");
    TomoFineGeometry.ReadNetCDF(tomogeometryfilename);

    //get the name of the file containing the data and read it in
    std::string tomodatafilename = jiba::AskFilename(
        "Tomography Data Filename: ");

    //read in data
    jiba::ReadTraveltimes(tomodatafilename, TomoData, TomoModel);
    //copy measurement configuration to refined model
    TomoFineGeometry.CopyMeasurementConfigurations(TomoModel);
    std::string scalgravdatafilename = jiba::AskFilename(
        "Scalar Gravity Data Filename: ");
    std::string ftgdatafilename = jiba::AskFilename("FTG Data Filename: ");
    std::string mtdatafilename = jiba::AskFilename("MT data filename: ");
    std::string gravmodelfilename = jiba::AskFilename(
        "Gravity Model Filename: ");
    std::string mtmodelfilename = jiba::AskFilename("MT Model Filename: ");
    jiba::ThreeDGravityModel GravModel;
    GravModel.ReadNetCDF(gravmodelfilename);
    //we only read in the gravity model for the information about the background layers
    //we take the actual gridding from the start model
    GravModel = TomoModel;

    jiba::ThreeDGravityModel::tMeasPosVec PosX, PosY, PosZ;
    jiba::ReadScalarGravityMeasurements(scalgravdatafilename, ScalGravData,
        PosX, PosY, PosZ);
    jiba::ReadTensorGravityMeasurements(ftgdatafilename, FTGData, PosX, PosY,
        PosZ);
    GravModel.ClearMeasurementPoints();
    for (size_t i = 0; i < PosX.size(); ++i)
      {
        GravModel.AddMeasurementPoint(PosX.at(i), PosY.at(i), PosZ.at(i));
      }

    //read in MT data
    std::vector<double> MTXPos, MTYPos, MTZPos, Frequencies;
    jiba::ReadImpedancesFromNetCDF(mtdatafilename, Frequencies, MTXPos, MTYPos,
        MTZPos, MTData);
    jiba::X3DModel MTModel;
    MTModel.ReadNetCDF(mtmodelfilename);
    //as for the gravity model the gridding is determined by the starting model
    //and we only read the mt model for the background layers
    MTModel = TomoModel;
    MTModel.ClearMeasurementPoints();
    std::copy(Frequencies.begin(), Frequencies.end(), std::back_inserter(
        MTModel.SetFrequencies()));
    for (size_t i = 0; i < MTXPos.size(); ++i)
      {
        MTModel.AddMeasurementPoint(MTXPos[i], MTYPos[i], MTZPos[i]);
      }

    //if we don't have data inversion doesn't make sense;
    if (TomoData.empty() || ScalGravData.empty() || MTData.empty())
      {
        std::cerr << "No measurements defined" << std::endl;
        exit(100);
      }

    jiba::rvec InvModel(TomoModel.GetSlownesses().num_elements());
    std::copy(TomoModel.GetSlownesses().origin(),
        TomoModel.GetSlownesses().origin()
            + TomoModel.GetSlownesses().num_elements(), InvModel.begin());

    const double z0 = 5.0;
    const double DepthExponent = -2.0;
    jiba::rvec WeightVector, ModelWeight(InvModel.size());
    //calculate the depth scaling
    jiba::ConstructDepthWeighting(GravModel.GetZCellSizes(), z0, WeightVector,
        jiba::WeightingTerm(DepthExponent));
    for (size_t i = 0; i < ModelWeight.size(); ++i)
      {
        ModelWeight(i) = WeightVector(i % GravModel.GetZCellSizes().size());
      }

    const double minslow = 1e-4;
    const double maxslow = 0.005;

    boost::shared_ptr<jiba::GeneralModelTransform> SlownessTransform(
        new jiba::TanhTransform(minslow, maxslow));
    boost::shared_ptr<jiba::GeneralModelTransform> DensityTransform(
        new jiba::DensityTransform(SlownessTransform));
    boost::shared_ptr<jiba::GeneralModelTransform> MTTransform(
        new jiba::ConductivityTransform(SlownessTransform));

    //double average = std::accumulate(InvModel.begin(),InvModel.end(),0.0)/InvModel.size();
    //std::fill(RefModel.begin(),RefModel.end(),average);
    jiba::rvec RefModel(InvModel);
    std::cout << "Background layers: "
        << GravModel.GetBackgroundDensities().size() << std::endl;


    boost::shared_ptr<jiba::TomographyObjective> TomoObjective(
        new jiba::TomographyObjective());
    TomoObjective->SetObservedData(TomoData);
    TomoObjective->SetFineModelGeometry(TomoFineGeometry);
    TomoObjective->SetCoarseModelGeometry(TomoModel);
    jiba::rvec TomoCovar(TomoData.size());
    //we assume a general error of 5 ms for the seismic data
    std::fill(TomoCovar.begin(), TomoCovar.end(), 5.0e-3);
    TomoObjective->SetDataCovar(TomoCovar);


    boost::shared_ptr<jiba::GravityObjective> ScalGravObjective(
        new jiba::GravityObjective(false, wantcuda));
    ScalGravObjective->SetObservedData(ScalGravData);
    ScalGravObjective->SetModelGeometry(GravModel);
    ScalGravObjective->SetDataCovar(jiba::ConstructError(ScalGravData, 0.0,
        5e-7));

    boost::shared_ptr<jiba::GravityObjective> FTGObjective(
        new jiba::GravityObjective(true, wantcuda));
    FTGObjective->SetObservedData(FTGData);
    FTGObjective->SetModelGeometry(GravModel);
    FTGObjective->SetDataCovar(jiba::ConstructError(FTGData, 0.02, 1e-9));


    boost::shared_ptr<jiba::X3DObjective> MTObjective(new jiba::X3DObjective());
    MTObjective->SetModelGeometry(MTModel);
    MTObjective->SetObservedData(MTData);
    MTObjective->SetDataCovar(jiba::ConstructMTError(MTData, 0.02));


    boost::shared_ptr<jiba::JointObjective> Objective(
        new jiba::JointObjective());

    boost::shared_ptr<jiba::MatOpRegularization> Regularization;
    if (wantcurvreg)
      {
        Regularization = boost::shared_ptr<jiba::MatOpRegularization>(
            new jiba::CurvatureRegularization(GravModel));
      }
    else
      {
        Regularization = boost::shared_ptr<jiba::MatOpRegularization>(
            new jiba::GradientRegularization(GravModel));
      }
    jiba::rvec ModCov(InvModel.size() * 3);
    //std::fill(ModCov.begin(), ModCov.end(), 1.0);
    ublas::vector_range<jiba::rvec>(ModCov, ublas::range(0, InvModel.size()))
        = InvModel;
    ublas::vector_range<jiba::rvec>(ModCov, ublas::range(InvModel.size(), 2
        * InvModel.size())) = InvModel;
    ublas::vector_range<jiba::rvec>(ModCov, ublas::range(2 * InvModel.size(), 3
        * InvModel.size())) = InvModel;

    if (vm.count("refmod"))
      {
        Regularization->SetReferenceModel(InvModel);
      }
    else
      {
        jiba::rvec ZeroMod(InvModel.size());
        ZeroMod.clear();
        Regularization->SetReferenceModel(ZeroMod);
      }
    //Regularization->SetDataCovar(ModCov);
    if (vm.count("xreg"))
      {
        const double xreg = vm["xreg"].as<double> ();
        std::cout << "Setting xreg: " << xreg << std::endl;
        Regularization->SetXWeight(xreg);
      }
    if (vm.count("yreg"))
      {
        const double yreg = vm["yreg"].as<double> ();
        std::cout << "Setting yreg: " << yreg << std::endl;
        Regularization->SetYWeight(yreg);
      }
    if (vm.count("zreg"))
      {
        const double zreg = vm["zreg"].as<double> ();
        std::cout << "Setting zreg: " << zreg << std::endl;
        Regularization->SetZWeight(zreg);
      }

    double tomolambda = 1.0;
    double scalgravlambda = 1.0;
    double ftglambda = 1.0;
    double reglambda = 1.0;
    double mtlambda = 1.0;
    std::cout << "Tomography Lambda: ";
    std::cin >> tomolambda;
    std::cout << "Scalar Gravimetry Lambda: ";
    std::cin >> scalgravlambda;
    std::cout << "FTG Lambda: ";
    std::cin >> ftglambda;
    std::cout << "MT Lambda: ";
    std::cin >> mtlambda;
    std::cout << "Regularization Lambda: ";
    std::cin >> reglambda;
    size_t maxiter = 1;
    std::cout << "Maximum iterations: ";
    std::cin >> maxiter;

    boost::posix_time::ptime starttime =
        boost::posix_time::microsec_clock::local_time();
    size_t ndata = 0;
    if (tomolambda > 0.0)
      {
        ndata += TomoData.size();
        Objective->AddObjective(TomoObjective, SlownessTransform, tomolambda);
        std::cout << "Tomo ndata: " << TomoData.size() << std::endl;
        std::cout << "Tomo lambda: " << tomolambda << std::endl;
      }
    if (scalgravlambda > 0.0)
      {
        ndata += ScalGravData.size();
        Objective->AddObjective(ScalGravObjective, DensityTransform,
            scalgravlambda);
        std::cout << "Scalar Gravity ndata: " << ScalGravData.size()
            << std::endl;
        std::cout << "Scalar Gravity lambda: " << scalgravlambda << std::endl;
      }
    if (ftglambda > 0.0)
      {
        ndata += FTGData.size();
        Objective->AddObjective(FTGObjective, DensityTransform, ftglambda);
        std::cout << "FTG ndata: " << FTGData.size() << std::endl;
        std::cout << "FTG lambda: " << ftglambda << std::endl;
      }
    if (mtlambda > 0.0)
      {
        ndata += MTData.size();
        Objective->AddObjective(MTObjective, MTTransform, mtlambda);
        std::cout << "MT ndata: " << MTData.size() << std::endl;
        std::cout << "MT lambda: " << mtlambda << std::endl;
      }
    Objective->AddObjective(Regularization, SlownessTransform, reglambda);
    std::cout << "Performing inversion." << std::endl;

    boost::shared_ptr<jiba::GradientBasedOptimization> Optimizer;
    if (wantnlcg)
      {
        Optimizer = boost::shared_ptr<jiba::GradientBasedOptimization>(
            new jiba::NonLinearConjugateGradient(Objective));
      }
    else
      {
        Optimizer = boost::shared_ptr<jiba::GradientBasedOptimization>(
            new jiba::LimitedMemoryQuasiNewton(Objective, correctionpairs));
      }

    //Optimizer->SetModelCovDiag(ModelWeight);

    size_t iteration = 0;
    InvModel = SlownessTransform->PhysicalToGeneralized(InvModel);
    jiba::rvec TomoInvModel(SlownessTransform->GeneralizedToPhysical(InvModel));
    std::ofstream misfitfile("misfit.out");
    //calculate initial misfit
    misfitfile << "0 " << Objective->CalcMisfit(InvModel) << " ";
    std::copy(Objective->GetIndividualFits().begin(),
        Objective->GetIndividualFits().end(), std::ostream_iterator<double>(
            misfitfile, " "));
    misfitfile << std::endl;
    bool terminate = true;
    do
      {
        terminate = true;
        try
          {
            std::cout << "Iteration: " << iteration << std::endl;
            std::copy(InvModel.begin(), InvModel.begin()
                + TomoModel.GetSlownesses().num_elements(),
                TomoModel.SetSlownesses().origin());
            TomoModel.WriteVTK("raw_model" + jiba::stringify(iteration)
                + ".tomo.inv.vtk");
            Optimizer->MakeStep(InvModel);

            ++iteration;

            TomoInvModel = SlownessTransform->GeneralizedToPhysical(InvModel);
            std::copy(TomoInvModel.begin(), TomoInvModel.begin()
                + TomoModel.GetSlownesses().num_elements(),
                TomoModel.SetSlownesses().origin());
            TomoModel.WriteVTK(modelfilename + jiba::stringify(iteration)
                + ".tomo.inv.vtk");
            std::cout << "Currrent Misfit: " << Optimizer->GetMisfit()
                << std::endl;
            std::cout << "Currrent Gradient: " << Optimizer->GetGradNorm()
                << std::endl;
            misfitfile << iteration << " " << Optimizer->GetMisfit() << " ";
            std::copy(Objective->GetIndividualFits().begin(),
                Objective->GetIndividualFits().end(), std::ostream_iterator<
                    double>(misfitfile, " "));
            std::cout << "Gradient Norms: ";
            std::copy(Objective->GetIndividualGradNorms().begin(),
                Objective->GetIndividualGradNorms().end(),
                std::ostream_iterator<double>(std::cout, " "));
            misfitfile << " " << Objective->GetNEval();
            misfitfile << std::endl;
            std::cout << "\n\n";
          } catch (jiba::FatalException &e)
          {
            std::cerr << e.what() << std::endl;
            iteration = maxiter;
          }

        for (size_t i = 0; i < Objective->GetIndividualFits().size() - 1; ++i)
          {
            if (Objective->GetIndividualFits().at(i) > Objective->GetObjective(
                i).GetNData())
              {
                terminate = false;
              }
            else
              {
                std::cout << "Reached target misfit." << std::endl;
                ;
                std::cout << "Objective number: " << i << std::endl;
                std::cout << "Misfit: " << Objective->GetIndividualFits().at(i)
                    << std::endl;
                std::cout << "Target: "
                    << Objective->GetObjective(i).GetNData() << std::endl;
              }
          }
        if (jiba::WantAbort())
          {
            terminate = true;
            jiba::RemoveAbort();
          }
      } while (iteration < maxiter && !terminate && Optimizer->GetGradNorm()
        > 1e-6);

    jiba::rvec DensInvModel(DensityTransform->GeneralizedToPhysical(InvModel));
    jiba::rvec CondInvModel(MTTransform->GeneralizedToPhysical(InvModel));
    std::copy(TomoInvModel.begin(), TomoInvModel.begin()
        + TomoModel.GetSlownesses().num_elements(),
        TomoModel.SetSlownesses().origin());
    std::copy(DensInvModel.begin(), DensInvModel.begin()
        + GravModel.SetDensities().num_elements(),
        GravModel.SetDensities().origin());
    std::copy(CondInvModel.begin(), CondInvModel.begin()
        + MTModel.GetConductivities().num_elements(),
        MTModel.SetConductivities().origin());

    //calculate the predicted refraction data
    std::cout << "Calculating response of inversion model." << std::endl;
    jiba::rvec TomoInvData(jiba::TomographyCalculator().Calculate(TomoModel));
    jiba::SaveTraveltimes(modelfilename + ".inv_tt.nc", TomoInvData, TomoModel);
    TomoModel.WriteNetCDF(modelfilename + ".tomo.inv.nc");
    TomoModel.WriteVTK(modelfilename + ".tomo.inv.vtk");

    //and write out the data and model
    //here we have to distinguish again between scalar and ftg data
    jiba::rvec ScalGravInvData(jiba::CreateGravityCalculator<
        jiba::MinMemGravityCalculator>::MakeScalar()->Calculate(GravModel));
    jiba::rvec FTGInvData(jiba::CreateGravityCalculator<
        jiba::MinMemGravityCalculator>::MakeTensor()->Calculate(GravModel));
    jiba::SaveScalarGravityMeasurements(modelfilename + ".inv_sgd.nc",
        ScalGravInvData, GravModel.GetMeasPosX(), GravModel.GetMeasPosY(),
        GravModel.GetMeasPosZ());
    jiba::SaveTensorGravityMeasurements(modelfilename + ".inv_ftg.nc",
        FTGInvData, GravModel.GetMeasPosX(), GravModel.GetMeasPosY(),
        GravModel.GetMeasPosZ());
    //calculate MT inversion result
    jiba::rvec MTInvData(jiba::X3DMTCalculator().Calculate(MTModel));
    jiba::WriteImpedancesToNetCDF(modelfilename + "inv_mt.nc",
        MTModel.GetFrequencies(), MTModel.GetMeasPosX(), MTModel.GetMeasPosY(),
        MTModel.GetMeasPosZ(), MTInvData);
    std::cout << "Writing out inversion results." << std::endl;

    GravModel.WriteVTK(modelfilename + ".grav.inv.vtk");
    GravModel.WriteNetCDF(modelfilename + ".grav.inv.nc");
    MTModel.WriteVTK(modelfilename + ".mt.inv.vtk");
    MTModel.WriteNetCDF(modelfilename + ".mt.inv.nc");
    std::ofstream datadiffile("data.diff");
    std::copy(Objective->GetDataDifference().begin(),
        Objective->GetDataDifference().end(), std::ostream_iterator<double>(
            datadiffile, "\n"));
    boost::posix_time::ptime endtime =
        boost::posix_time::microsec_clock::local_time();
    double cachedruntime = (endtime - starttime).total_seconds();
    std::cout << "Runtime: " << cachedruntime << " s" << std::endl;
    std::cout << std::endl;
  }