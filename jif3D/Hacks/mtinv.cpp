//============================================================================
// Name        : gravinv.cpp
// Author      : Max Moorkamp
// Version     :
// Copyright   : 2008, MM
//============================================================================

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <cstdlib>
#ifdef HAVEOPENMP
#include <omp.h>
#endif

#include <boost/program_options.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/log/core.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/expressions.hpp>

#include "../Global/convert.h"
#include "../Global/FatalException.h"
#include "../Global/NumUtil.h"
#include "../Global/VectorTransform.h"
#include "../Global/FileUtil.h"
#include "../Global/Noise.h"
#include "../Global/ReadWriteSparseMatrix.h"
#include "../ModelBase/VTKTools.h"
#include "../ModelBase/NetCDFModelTools.h"
#include "../ModelBase/ReadAnyModel.h"
#include "../Regularization/GradientRegularization.h"
#include "../Regularization/CrossGradient.h"
#include "../Inversion/LimitedMemoryQuasiNewton.h"
#include "../Inversion/JointObjective.h"
#include "../Inversion/ModelTransforms.h"
#include "../Inversion/ThreeDModelObjective.h"
#include "../Regularization/MinDiffRegularization.h"
#include "../MT/X3DModel.h"
#include "../MT/X3DMTCalculator.h"
#include "../MT/ReadWriteImpedances.h"
#include "../Joint/SetupRegularization.h"
#include "../Joint/SetupInversion.h"
#include "../Joint/SetupMT.h"
#include "../Joint/InversionOutput.h"

namespace ublas = boost::numeric::ublas;
namespace po = boost::program_options;
namespace logging = boost::log;

int main(int argc, char *argv[])
  {

    double mincond = 1e-6;
    double maxcond = 10;
    double relerr = 0.02;
    double coolingfactor = 1.0;
    double xorigin = 0.0;
    double yorigin = 0.0;
    std::string X3DName = "x3d";
    std::string MTInvCovarName;
    std::string RefModelName;
    std::string CrossModelName;
    boost::shared_ptr<jif3D::JointObjective> Objective(new jif3D::JointObjective(true));

    double DistCorr = 0;
    po::options_description desc("General options");
    desc.add_options()("help", "produce help message")("threads", po::value<int>(),
        "The number of openmp threads")("xorigin",
        po::value(&xorigin)->default_value(0.0),
        "The origin for the inversion grid in x-direction")("yorigin",
        po::value(&yorigin)->default_value(0.0),
        "The origin for the inversion grid in y-direction")("covmod",
        po::value<std::string>(), "A file containing the model covariance")("mincond",
        po::value(&mincond)->default_value(1e-4),
        "The minimum value for conductivity in S/m")("maxcond",
        po::value(&maxcond)->default_value(5),
        "The maximum value for conductivity in S/m")("tempdir", po::value<std::string>(),
        "The name of the directory to store temporary files in")("distcorr",
        po::value(&DistCorr)->default_value(0), "Correct for distortion within inversion")(
        "x3dname", po::value(&X3DName)->default_value("x3d"),
        "The name of the executable for x3d")("mtinvcovar",
        po::value<std::string>(&MTInvCovarName),
        "Inverse covariance matrix to use in MT misfit calculation.")("inderrors",
        "Use the individual errors for each element instead of the same for all elements")(
        "mtrelerr", po::value(&relerr)->default_value(0.02),
        "Error floor for impedance estimates.")("coolingfactor",
        po::value(&coolingfactor)->default_value(1.0),
        "The factor to multiply the weight for the regularization at each iteration EXPERIMENTAL")(
        "debug", "Show debugging output.")("opt",
        "Use opt for Green's function calculation in x3d.")("refmodel",
        po::value(&RefModelName),
        "The name of the reference model to substract before calculating smoothness")(
        "crossmodel", po::value(&CrossModelName),
        "The name of a model to use as a cross-gradient constraint");

    jif3D::SetupRegularization RegSetup;
    jif3D::SetupInversion InversionSetup;

    desc.add(RegSetup.SetupOptions());
    desc.add(InversionSetup.SetupOptions());
    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);

    if (vm.count("debug"))
      {
        logging::core::get()->set_filter(
            logging::trivial::severity >= logging::trivial::debug);
      }
    else
      {
        logging::core::get()->set_filter(
            logging::trivial::severity >= logging::trivial::warning);
      }
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
    jif3D::rvec CovModVec;
    if (vm.count("covmod"))
      {
        jif3D::ThreeDModelBase CovModel;
        //we store the covariances in a seismic model file
        //but we do not have to have an equidistant grid
        std::string Filename(vm["covmod"].as<std::string>());
        CovModel = *jif3D::ReadAnyModel(Filename).get();
        const size_t ncovmod = CovModel.GetData().num_elements();
        CovModVec.resize(ncovmod);
        std::copy(CovModel.GetData().origin(), CovModel.GetData().origin() + ncovmod,
            CovModVec.begin());
        BOOST_LOG_TRIVIAL(debug)<< "Reading in covariance file: " << Filename << std::endl;
        BOOST_LOG_TRIVIAL(debug)<< "Covariance file has: " << CovModel.GetNModelElements() << " elements" << std::endl;
        BOOST_LOG_TRIVIAL(debug)<< "Covariance file has: " << std::count(CovModel.GetData().origin(),CovModel.GetData().origin() + ncovmod,1) << " elements = 1" << std::endl;
      }

    boost::filesystem::path TempDir = boost::filesystem::current_path();
    if (vm.count("tempdir"))
      {
        TempDir = vm["tempdir"].as<std::string>();
        if (!boost::filesystem::is_directory(TempDir))
          {
            std::cerr << TempDir.string() << " is not a directory or does not exist ! \n";
            return 500;
          }
      }

//first we read in the starting model and the measured data
    std::string modelfilename = jif3D::AskFilename("Starting model Filename: ");
    std::string extension = jif3D::GetFileExtension(modelfilename);
    std::cout << "Extension: " << extension << std::endl;
//we read in the starting modelfile
    jif3D::X3DModel Model;
    Model.SetOrigin(xorigin, yorigin, 0.0);
    if (extension.compare(".mod") == 0)
      {
        Model.ReadModEM(modelfilename);
        Model.WriteNetCDF(modelfilename + ".nc");
      }
    else
      {
        Model.ReadNetCDF(modelfilename);
        Model.WriteModEM(modelfilename + ".dat");
      }
//get the name of the file containing the data and read it in
    std::string datafilename = jif3D::AskFilename("Data Filename: ");
    extension = jif3D::GetFileExtension(datafilename);
//read in data
    jif3D::rvec Data, ZError;
    std::vector<double> XCoord, YCoord, ZCoord, Frequencies, C;
    if (extension.compare(".dat") == 0)
      {
        jif3D::ReadImpedancesFromModEM(datafilename, Frequencies, XCoord, YCoord, ZCoord,
            Data, ZError);
        jif3D::WriteImpedancesToNetCDF("start_data.nc", Frequencies, XCoord, YCoord,
            ZCoord, Data, ZError);
      }
    else
      {
        jif3D::ReadImpedancesFromNetCDF(datafilename, Frequencies, XCoord, YCoord, ZCoord,
            Data, ZError, C);
        jif3D::WriteImpedancesToModEM("start_data.dat", Frequencies, XCoord, YCoord,
            ZCoord, Data, ZError);
      }
    //if we have not specified a background in the model file
    if (Model.GetBackgroundConductivities().empty())
      {
        jif3D::rvec Cond(Model.GetConductivities().num_elements());
        std::copy(Model.GetConductivities().origin(),
            Model.GetConductivities().origin() + Model.GetConductivities().num_elements(),
            Cond.begin());
        std::nth_element(Cond.begin(), Cond.begin() + Cond.size() / 2, Cond.end());
        double MedCond = Cond(Cond.size() / 2);
        std::vector<double> bg_thick(2), bg_cond(2, MedCond);
        bg_thick.at(0) = Model.GetZCoordinates().at(Model.GetZCoordinates().size() - 1);
        bg_thick.at(1) = bg_thick.at(0);
        Model.SetBackgroundConductivities(bg_cond);
        Model.SetBackgroundThicknesses(bg_thick);
      }
    std::copy(Frequencies.begin(), Frequencies.end(),
        std::back_inserter(Model.SetFrequencies()));
//if we don't have data inversion doesn't make sense;
    if (Data.empty())
      {
        std::cerr << "No measurements defined" << std::endl;
        exit(100);
      }
    for (size_t i = 0; i < XCoord.size(); ++i)
      {
        Model.AddMeasurementPoint(XCoord[i], YCoord[i], ZCoord[i]);
      }

    jif3D::rvec FirstFreq(XCoord.size());
    std::iota(FirstFreq.begin(), FirstFreq.end(), 1);
    //std::copy(Data.begin(), Data.begin() + XCoord.size(), FirstFreq.begin());
    jif3D::Write3DDataToVTK(datafilename + ".vtk", "MTSites", FirstFreq, XCoord, YCoord,
        ZCoord);

    Model.SetDistortionParameters(C);
//we define a few constants that are used throughout the inversion
    const size_t ndata = Data.size();

    for (size_t i = 0; i < Model.GetConductivities().shape()[2]; ++i)
      {
        Model.SetConductivities()[0][0][i] *= (1 + 0.01 * (i + 1));
      }

    const size_t ngrid = Model.GetConductivities().num_elements();
    jif3D::rvec InvModel(ngrid);
    std::copy(Model.GetConductivities().origin(),
        Model.GetConductivities().origin() + ngrid, InvModel.begin());
    Model.WriteVTK("start.vtk");

    //if we want to correct for distortion within the inversion
    if (DistCorr > 0)
      {
        //if we did not read distortion values together with the
        //observed data, we set the distortion matrix C at each site
        // to identity matrix
        if (C.empty())
          {
            C.resize(XCoord.size() * 4);
            for (size_t i = 0; i < XCoord.size(); ++i)
              {
                C[i * 4] = 1.0;
                C[i * 4 + 1] = 0.0;
                C[i * 4 + 2] = 0.0;
                C[i * 4 + 3] = 1.0;
              }
          }
        //we need to expand the model vector to hold the
        //elements of the distortion matrix
        jif3D::rvec Grid(InvModel);
        InvModel.resize(ngrid + C.size());
        std::copy(Grid.begin(), Grid.end(), InvModel.begin());
        std::copy(C.begin(), C.end(), InvModel.begin() + ngrid);
        //also the diagonal of the model covariance needs to
        //accommodate the new parameters
        jif3D::rvec OldCov(CovModVec);
        CovModVec.resize(InvModel.size());
        std::fill(CovModVec.begin(), CovModVec.end(), 1.0);
        std::copy(OldCov.begin(), OldCov.end(), CovModVec.begin());
      }

    boost::shared_ptr<jif3D::GeneralModelTransform> Copier(new jif3D::ModelCopyTransform);
    if (vm.count("crossmodel"))
      {
        boost::shared_ptr<jif3D::ThreeDModelBase> CrossModel(
            jif3D::ReadAnyModel(CrossModelName));
        if (CrossModel->GetNModelElements() != ngrid)
          {
            std::cerr
                << "Cross gradient model does not have the same size as inversion grid"
                << std::endl;
            return 100;
          }
        double cgweight = 1.0;
        std::cout << "Cross-gradient weight: ";
        std::cin >> cgweight;

        jif3D::rvec OldInv(InvModel);
        InvModel.resize(OldInv.size() + ngrid, true);
        CovModVec.resize(InvModel.size(), true);
        std::copy(CrossModel->GetData().origin(), CrossModel->GetData().origin() + ngrid,
            InvModel.begin() + OldInv.size());
        std::fill_n(CovModVec.begin(), OldInv.size(), 1.0);
        std::fill_n(CovModVec.begin() + OldInv.size(), ngrid, 1e-32);

        std::cout << "Inversion vector length: " << InvModel.size() << " " << ngrid << " "
            << OldInv.size() << std::endl;
        boost::shared_ptr<jif3D::MultiSectionTransform> CrossTrans(
            new jif3D::MultiSectionTransform(InvModel.size(), 0, ngrid, Copier));
        CrossTrans->AddSection(OldInv.size(), InvModel.size(), Copier);

        boost::shared_ptr<jif3D::CrossGradient> CrossConstr(
            new jif3D::CrossGradient(Model));
        Objective->AddObjective(CrossConstr, CrossTrans, cgweight, "Cross",
            jif3D::JointObjective::coupling);
      }

    boost::shared_ptr<jif3D::ChainedTransform> ConductivityTransform(
        new jif3D::ChainedTransform);
    jif3D::rvec RefModel(ngrid, 1.0);
//because the tanh transform is used inside a logarithmic transform
//we need to take the natural logarithm of the actual minimum and maximum
    ConductivityTransform->AppendTransform(
        boost::shared_ptr<jif3D::GeneralModelTransform>(
            new jif3D::TanhTransform(std::log(mincond), std::log(maxcond))));
    ConductivityTransform->AppendTransform(
        boost::shared_ptr<jif3D::GeneralModelTransform>(
            new jif3D::LogTransform(RefModel)));

    boost::shared_ptr<jif3D::MultiSectionTransform> MTTransform(
        new jif3D::MultiSectionTransform(InvModel.size(), 0, ngrid,
            ConductivityTransform));
    if (DistCorr > 0)
      {
        MTTransform->AddSection(ngrid, ngrid + C.size(), Copier);
      }

    subrange(InvModel, 0, ngrid) = subrange(MTTransform->PhysicalToGeneralized(InvModel),
        0, ngrid);

    bool WantDistCorr = (DistCorr > 0);
    jif3D::X3DMTCalculator Calculator(TempDir, X3DName, WantDistCorr);
    if (vm.count("opt"))
      {
        BOOST_LOG_TRIVIAL(info)<< "Using Opt type Green's functions ";
        Calculator.SetGreenType1(jif3D::GreenCalcType::opt);
        Calculator.SetGreenType4(jif3D::GreenCalcType::opt);
      }
    boost::shared_ptr<jif3D::ThreeDModelObjective<jif3D::X3DMTCalculator> > X3DObjective(
        new jif3D::ThreeDModelObjective<jif3D::X3DMTCalculator>(Calculator));

    X3DObjective->SetObservedData(Data);
    X3DObjective->SetCoarseModelGeometry(Model);

//create objects for the misfit and error estimates
    if (vm.count("mtinvcovar"))
      {
        jif3D::comp_mat InvCov(Data.size(), Data.size());
        jif3D::ReadSparseMatrixFromNetcdf(MTInvCovarName, InvCov, "InvCovariance");
        X3DObjective->SetInvCovMat(InvCov);
      }
    else
      {
        jif3D::rvec MinErr(ZError.size());
        if (vm.count("inderrors"))
          {
            MinErr = jif3D::ConstructError(Data, ZError, relerr);
          }
        else
          {
            MinErr = jif3D::ConstructMTError(Data, relerr);
          }
        for (size_t i = 0; i < ZError.size(); ++i)
          {
            ZError(i) = std::max(ZError(i), MinErr(i));

          }
        X3DObjective->SetDataError(ZError);
      }

    jif3D::rvec GridCov(ngrid, 1.0);
    //std::copy(CovModVec.begin(), CovModVec.begin() + ngrid, GridCov.begin());
    boost::shared_ptr<jif3D::RegularizationFunction> Regularization =
        RegSetup.SetupObjective(vm, Model, GridCov);

    double lambda = 1.0;
    std::cout << "Lambda: ";
    std::cin >> lambda;
    Objective->AddObjective(X3DObjective, MTTransform, 1.0, "MT",
        jif3D::JointObjective::datafit);
    boost::shared_ptr<jif3D::MultiSectionTransform> ModRegTrans(
        new jif3D::MultiSectionTransform(InvModel.size(), 0, ngrid, Copier));
    Objective->AddObjective(Regularization, ModRegTrans, lambda, "Regularization",
        jif3D::JointObjective::regularization);
    if (vm.count("refmodel"))
      {
        jif3D::X3DModel RefModel;
        RefModel.ReadNetCDF(RefModelName);
        jif3D::rvec RefVec(RefModel.GetNModelElements());
        std::copy(RefModel.GetConductivities().origin(),
            RefModel.GetConductivities().origin() + RefModel.GetNModelElements(),
            RefVec.begin());
        Regularization->SetReferenceModel(
            ConductivityTransform->PhysicalToGeneralized(RefVec));
      }

    if (DistCorr > 0)
      {
        jif3D::rvec CRef(XCoord.size() * 4);

        for (size_t i = 0; i < XCoord.size(); ++i)
          {
            CRef(i * 4) = 1.0;
            CRef(i * 4 + 1) = 0.0;
            CRef(i * 4 + 2) = 0.0;
            CRef(i * 4 + 3) = 1.0;
          }
        jif3D::X3DModel DistModel;
        DistModel.SetMeshSize(XCoord.size() * 4, 1, 1);
        boost::shared_ptr<jif3D::RegularizationFunction> DistReg(
            new jif3D::MinDiffRegularization(DistModel));
        DistReg->SetReferenceModel(CRef);
        boost::shared_ptr<jif3D::MultiSectionTransform> DistRegTrans(
            new jif3D::MultiSectionTransform(InvModel.size(), ngrid, ngrid + CRef.size(),
                Copier));

        Objective->AddObjective(DistReg, DistRegTrans, DistCorr, "DistReg",
            jif3D::JointObjective::regularization);
      }

    size_t maxiter = 30;
    std::cout << "Maximum number of iterations: ";
    std::cin >> maxiter;
    std::cout << "Performing inversion." << std::endl;

    boost::shared_ptr<jif3D::GradientBasedOptimization> Optimizer =
        InversionSetup.ConfigureInversion(vm, Objective, InvModel, CovModVec);

    std::ofstream misfitfile("misfit.out");
    std::ofstream rmsfile("rms.out");
    std::ofstream weightfile("weights.out");
    double InitialMisfit = Objective->CalcMisfit(InvModel);

    StoreMisfit(misfitfile, 0, InitialMisfit, *Objective);
    StoreRMS(rmsfile, 0, *Objective);

    size_t iteration = 0;
    boost::posix_time::ptime starttime = boost::posix_time::microsec_clock::local_time();

    bool terminate = false;
    while (iteration < maxiter && !terminate)
      {
        try
          {

            std::cout << "\n\n Iteration: " << iteration << std::endl;

            //update the inversion model
            Optimizer->MakeStep(InvModel);
            Objective->MultiplyWeights(jif3D::JointObjective::regularization,
                coolingfactor);
            ++iteration;

            //write out some information about misfit to the screen
            jif3D::rvec DistModel = MTTransform->GeneralizedToPhysical(InvModel);
            std::copy(DistModel.begin(), DistModel.begin() + Model.GetNModelElements(),
                Model.SetData().origin());
            Model.WriteVTK(modelfilename + jif3D::stringify(iteration) + ".mt.inv.vtk");
            Model.WriteNetCDF(modelfilename + jif3D::stringify(iteration) + ".mt.inv.nc");
            if (DistCorr > 0)
              {
                std::copy(DistModel.begin() + Model.GetNModelElements(),
                    DistModel.begin() + Model.GetNModelElements() + C.size(), C.begin());
                jif3D::WriteImpedancesToNetCDF(
                    modelfilename + jif3D::stringify(iteration) + ".dist_imp.nc",
                    Frequencies, XCoord, YCoord, ZCoord, Data, ZError, C);
              }
            std::cout << "Currrent Misfit: " << Optimizer->GetMisfit() << std::endl;
            std::cout << "Currrent Gradient: " << Optimizer->GetGradNorm() << std::endl;
            //and write the current misfit for all objectives to a misfit file
            StoreMisfit(misfitfile, iteration, Optimizer->GetMisfit(), *Objective);
            StoreRMS(rmsfile, iteration, *Objective);
            StoreWeights(weightfile, iteration, *Objective);
            std::cout << "\n\n";

            terminate = jif3D::WantAbort();
          } catch (jif3D::FatalException &e)
          {
            std::cerr << e.what() << std::endl;
            terminate = true;
          }
      }

    boost::posix_time::ptime endtime = boost::posix_time::microsec_clock::local_time();
    double cachedruntime = (endtime - starttime).total_seconds();
    std::cout << "Runtime: " << cachedruntime << " s" << std::endl;
    std::cout << std::endl;

    if (vm.count("crossmodel"))
      {
        std::copy(InvModel.end() - Model.GetNModelElements(), InvModel.end(),
            Model.SetConductivities().origin());
        Model.WriteVTK("crossref.final.vtk");
      }

    InvModel = MTTransform->GeneralizedToPhysical(InvModel);

    std::copy(InvModel.begin(),
        InvModel.begin() + Model.GetConductivities().num_elements(),
        Model.SetConductivities().origin());
    std::copy(InvModel.begin() + Model.GetNModelElements(), InvModel.end(), C.begin());

    Model.SetDistortionParameters(C);
//calculate the predicted data
    std::cout << "Calculating response of inversion model." << std::endl;
    jif3D::rvec InvData(X3DObjective->GetSyntheticData());
    jif3D::WriteImpedancesToNetCDF(modelfilename + ".inv_imp.nc", Frequencies, XCoord,
        YCoord, ZCoord, InvData, X3DObjective->GetDataError(), C);
    jif3D::WriteImpedancesToNetCDF(modelfilename + ".dist_imp.nc", Frequencies, XCoord,
        YCoord, ZCoord, Data, ZError, C);
    jif3D::WriteImpedancesToNetCDF(modelfilename + ".diff_imp.nc", Frequencies, XCoord,
        YCoord, ZCoord, X3DObjective->GetIndividualMisfit());

//and write out the data and model
//here we have to distinguish again between scalar and ftg data
    std::cout << "Writing out inversion results." << std::endl;

    Model.WriteVTK(modelfilename + ".inv.vtk");
    Model.WriteNetCDF(modelfilename + ".inv.nc");



    std::cout << std::endl;
  }
