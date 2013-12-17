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

#include "../Global/convert.h"
#include "../Global/FatalException.h"
#include "../Global/NumUtil.h"
#include "../Global/VectorTransform.h"
#include "../Global/FileUtil.h"
#include "../Global/Noise.h"
#include "../ModelBase/VTKTools.h"
#include "../ModelBase/NetCDFModelTools.h"
#include "../Regularization/GradientRegularization.h"
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

int main(int argc, char *argv[])
  {

    double mincond = 1e-6;
    double maxcond = 10;

    boost::shared_ptr<jif3D::JointObjective> Objective(new jif3D::JointObjective(true));

    bool WantDistCorr;
    po::options_description desc("General options");
    desc.add_options()("help", "produce help message")("threads", po::value<int>(),
        "The number of openmp threads")("covmod", po::value<std::string>(),
        "A file containing the model covariance")("mincond",
        po::value(&mincond)->default_value(1e-4),
        "The minimum value for conductivity in S/m")("maxcond",
        po::value(&maxcond)->default_value(5),
        "The maximum value for conductivity in S/m")("tempdir", po::value<std::string>(),
        "The name of the directory to store temporary files in")("distcorr",
        po::value(&WantDistCorr)->default_value(false),
        "Correct for distortion within inversion");

    jif3D::SetupRegularization RegSetup;
    jif3D::SetupInversion InversionSetup;
    jif3D::SetupMT MTSetup;
    desc.add(MTSetup.SetupOptions());
    desc.add(RegSetup.SetupOptions());
    desc.add(InversionSetup.SetupOptions());
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
    jif3D::rvec CovModVec;
    if (vm.count("covmod"))
      {
        jif3D::X3DModel CovModel;
        CovModel.ReadNetCDF(vm["covmod"].as<std::string>());
        const size_t ncovmod = CovModel.GetConductivities().num_elements();
        CovModVec.resize(ncovmod);
        std::copy(CovModel.GetConductivities().origin(),
            CovModel.GetConductivities().origin() + ncovmod, CovModVec.begin());
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

    //we read in the starting modelfile
    jif3D::X3DModel Model;
    Model.ReadNetCDF(modelfilename);

    //get the name of the file containing the data and read it in
    std::string datafilename = jif3D::AskFilename("Data Filename: ");

    //read in data
    jif3D::rvec Data, ZError;
    std::vector<double> XCoord, YCoord, ZCoord, Frequencies, C;
    jif3D::ReadImpedancesFromNetCDF(datafilename, Frequencies, XCoord, YCoord, ZCoord,
        Data, ZError, C);
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
    std::copy(Data.begin(), Data.begin() + XCoord.size(), FirstFreq.begin());
    jif3D::Write3DDataToVTK(datafilename + ".vtk", "Z", FirstFreq, XCoord, YCoord,
        ZCoord);

    Model.SetDistortionParameters(C);
    //we define a few constants that are used throughout the inversion
    const size_t ndata = Data.size();
    const double errorlevel = 0.02;

    //create objects for the misfit and a very basic error estimate
    jif3D::rvec MTDataError = jif3D::ConstructMTError(Data, errorlevel);
    std::transform(ZError.begin(), ZError.end(), MTDataError.begin(), ZError.begin(),
        [] (double a, double b)
          { return std::max(a,b);});

    for (size_t i = 0; i < Model.GetConductivities().shape()[2]; ++i)
      {
        Model.SetConductivities()[0][0][i] *= (1 + 0.0001 * (i + 1));
      }

    const size_t ngrid = Model.GetConductivities().num_elements();
    jif3D::rvec InvModel(ngrid);
    std::copy(Model.GetConductivities().origin(),
        Model.GetConductivities().origin() + ngrid, InvModel.begin());
    Model.WriteVTK("start.vtk");

    if (WantDistCorr)
      {
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
        jif3D::rvec Grid(InvModel);
        InvModel.resize(ngrid + C.size());
        std::copy(Grid.begin(), Grid.end(), InvModel.begin());
        std::copy(C.begin(), C.end(), InvModel.begin() + ngrid);
      }

    boost::shared_ptr<jif3D::ChainedTransform> ConductivityTransform(
        new jif3D::ChainedTransform);
    jif3D::rvec RefModel(InvModel);
    std::fill(RefModel.begin(), RefModel.end(), 1.0);
    //because the tanh transform is used inside a logarithmic transform
    //we need to take the natural logarithm of the actual minimum and maximum
    ConductivityTransform->AppendTransform(
        boost::shared_ptr<jif3D::GeneralModelTransform>(
            new jif3D::TanhTransform(std::log(mincond), std::log(maxcond))));
    ConductivityTransform->AppendTransform(
        boost::shared_ptr<jif3D::GeneralModelTransform>(
            new jif3D::LogTransform(RefModel)));

    boost::shared_ptr<jif3D::GeneralModelTransform> Copier(new jif3D::ModelCopyTransform);
    boost::shared_ptr<jif3D::MultiSectionTransform> MTTransform(
        new jif3D::MultiSectionTransform(InvModel.size(), 0, ngrid,
            ConductivityTransform));
    if (WantDistCorr)
      {
        MTTransform->AddSection(ngrid, InvModel.size(), Copier);
      }

    InvModel = MTTransform->PhysicalToGeneralized(InvModel);

    jif3D::X3DMTCalculator Calculator(TempDir, WantDistCorr);
    boost::shared_ptr<jif3D::ThreeDModelObjective<jif3D::X3DMTCalculator> > X3DObjective(
        new jif3D::ThreeDModelObjective<jif3D::X3DMTCalculator>(Calculator));

    X3DObjective->SetObservedData(Data);
    X3DObjective->SetCoarseModelGeometry(Model);
    X3DObjective->SetDataError(ZError);

    boost::shared_ptr<jif3D::RegularizationFunction> Regularization =
        RegSetup.SetupObjective(vm, Model, CovModVec);

    double lambda = 1.0;
    std::cout << "Lambda: ";
    std::cin >> lambda;
    Objective->AddObjective(X3DObjective, MTTransform, 1.0, "MT",
        jif3D::JointObjective::datafit);
    Objective->AddObjective(Regularization, Copier, lambda, "Regularization",
        jif3D::JointObjective::regularization);

    if (WantDistCorr)
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
            new jif3D::MultiSectionTransform(InvModel.size(), ngrid, InvModel.size(),
                Copier));
        double distreglambda = 0;

        std::cout << "Regularization for distortion: ";
        std::cin >> distreglambda;
        Objective->AddObjective(DistReg, DistRegTrans, distreglambda, "DistReg",
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
    double InitialMisfit = Objective->CalcMisfit(InvModel);

    StoreMisfit(misfitfile, 0, InitialMisfit, *Objective);
    StoreRMS(rmsfile, 0, *Objective);

    size_t iteration = 0;
    boost::posix_time::ptime starttime = boost::posix_time::microsec_clock::local_time();

    bool terminate = false;
    while (iteration < maxiter && !terminate && Objective->GetIndividualFits()[0] > ndata)
      {
        try
          {

            std::cout << "\n\n Iteration: " << iteration << std::endl;

            //update the inversion model
            Optimizer->MakeStep(InvModel);

            ++iteration;

            SaveModel(InvModel, *ConductivityTransform.get(), Model,
                modelfilename + jif3D::stringify(iteration) + ".mt.inv");
            //write out some information about misfit to the screen
            std::cout << "Currrent Misfit: " << Optimizer->GetMisfit() << std::endl;
            std::cout << "Currrent Gradient: " << Optimizer->GetGradNorm() << std::endl;
            //and write the current misfit for all objectives to a misfit file
            StoreMisfit(misfitfile, iteration, Optimizer->GetMisfit(), *Objective);
            StoreRMS(rmsfile, iteration, *Objective);
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
    InvModel = ConductivityTransform->GeneralizedToPhysical(InvModel);

    std::copy(InvModel.begin(),
        InvModel.begin() + Model.GetConductivities().num_elements(),
        Model.SetConductivities().origin());
    std::copy(InvModel.begin() + Model.GetNModelElements(), InvModel.end(), C.begin());
    Model.SetDistortionParameters(C);
    //calculate the predicted data
    std::cout << "Calculating response of inversion model." << std::endl;
    jif3D::rvec InvData(jif3D::X3DMTCalculator().Calculate(Model));
    jif3D::WriteImpedancesToNetCDF(modelfilename + ".inv_imp.nc", Frequencies, XCoord,
        YCoord, ZCoord, InvData, X3DObjective->GetDataError(), C);
    jif3D::WriteImpedancesToNetCDF(modelfilename + ".diff_imp.nc", Frequencies, XCoord,
        YCoord, ZCoord, X3DObjective->GetIndividualMisfit());

    //and write out the data and model
    //here we have to distinguish again between scalar and ftg data
    std::cout << "Writing out inversion results." << std::endl;

    Model.WriteVTK(modelfilename + ".inv.vtk");
    Model.WriteNetCDF(modelfilename + ".inv.nc");
    std::cout << std::endl;
  }
