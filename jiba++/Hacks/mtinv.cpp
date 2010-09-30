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
#include <omp.h>
#include <boost/bind.hpp>
#include <boost/program_options.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include "../Global/convert.h"
#include "../Global/FatalException.h"
#include "../Global/NumUtil.h"
#include "../Global/VectorTransform.h"
#include "../Global/FileUtil.h"
#include "../Global/Noise.h"
#include "../ModelBase/VTKTools.h"
#include "../ModelBase/NetCDFTools.h"
#include "../Regularization/GradientRegularization.h"
#include "../Inversion/LimitedMemoryQuasiNewton.h"
#include "../Inversion/JointObjective.h"
#include "../Inversion/ModelTransforms.h"
#include "../MT/X3DModel.h"
#include "../MT/X3DObjective.h"
#include "../MT/X3DMTCalculator.h"
#include "../MT/ReadWriteImpedances.h"
#include "../Joint/SetupRegularization.h"
#include "../Joint/SetupInversion.h"

namespace ublas = boost::numeric::ublas;
namespace po = boost::program_options;

int main(int argc, char *argv[])
  {

    double mincond = 1e-6;
    double maxcond = 10;

    po::options_description desc("General options");
    desc.add_options()("help", "produce help message")("threads",
        po::value<int>(), "The number of openmp threads")("covmod", po::value<
        std::string>(), "A file containing the model covariance")("mincond",
        po::value(&mincond)->default_value(1e-6),
        "The minimum value for conductivity in S/m")("maxcond", po::value(
        &maxcond)->default_value(10),
        "The maximum value for conductivity in S/m");

    jiba::SetupRegularization RegSetup;
    jiba::SetupInversion InversionSetup;
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

    if (vm.count("threads"))
      {
        omp_set_num_threads(vm["threads"].as<int> ());
      }

    jiba::rvec CovModVec;
    if (vm.count("covmod"))
      {
        jiba::X3DModel CovModel;
        CovModel.ReadNetCDF(vm["covmod"].as<std::string> ());
        const size_t ncovmod = CovModel.GetConductivities().num_elements();
        CovModVec.resize(ncovmod);
        std::copy(CovModel.GetConductivities().origin(),
            CovModel.GetConductivities().origin() + ncovmod, CovModVec.begin());
      }

    //first we read in the starting model and the measured data
    std::string modelfilename = jiba::AskFilename("Starting model Filename: ");

    //we read in the starting modelfile
    jiba::X3DModel Model;
    Model.ReadNetCDF(modelfilename);

    //get the name of the file containing the data and read it in
    std::string datafilename = jiba::AskFilename("Data Filename: ");

    //read in data
    jiba::rvec Data;
    std::vector<double> XCoord, YCoord, ZCoord, Frequencies;
    jiba::ReadImpedancesFromNetCDF(datafilename, Frequencies, XCoord, YCoord,
        ZCoord, Data);
    std::copy(Frequencies.begin(), Frequencies.end(), std::back_inserter(
        Model.SetFrequencies()));
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

    jiba::rvec FirstFreq(XCoord.size());
    std::copy(Data.begin(), Data.begin() + XCoord.size(), FirstFreq.begin());
    jiba::Write3DDataToVTK(datafilename + ".vtk", "Z", FirstFreq, XCoord,
        YCoord, ZCoord);

    //we define a few constants that are used throughout the inversion
    const size_t ndata = Data.size();
    const double errorlevel = 0.02;

    //create objects for the misfit and a very basic error estimate
    jiba::rvec MTDataError = jiba::ConstructMTError(Data, errorlevel);
    std::ofstream errorfile("error.out");
    std::copy(MTDataError.begin(), MTDataError.end(), std::ostream_iterator<
        double>(errorfile, "\n"));
    std::flush(errorfile);

    for (size_t i = 0; i < Model.GetConductivities().shape()[2]; ++i)
      {
        Model.SetConductivities()[0][0][i] *= (1 + 0.0001 * (i + 1));
      }
    jiba::rvec InvModel(Model.GetConductivities().num_elements());
    std::copy(Model.GetConductivities().origin(),
        Model.GetConductivities().origin()
            + Model.GetConductivities().num_elements(), InvModel.begin());
    Model.WriteVTK("start.vtk");

    boost::shared_ptr<jiba::ChainedTransform> ConductivityTransform(
        new jiba::ChainedTransform);

    jiba::rvec RefModel(InvModel);
    std::fill(RefModel.begin(), RefModel.end(), 1.0);
    //because the tanh transform is used inside a logarithmic transform
    //we need to take the natural logarithm of the actual minimum and maximum
    ConductivityTransform->AddTransform(boost::shared_ptr<
        jiba::GeneralModelTransform>(new jiba::TanhTransform(std::log(mincond),
        std::log(maxcond))));
    ConductivityTransform->AddTransform(boost::shared_ptr<
        jiba::GeneralModelTransform>(new jiba::LogTransform(RefModel)));

    InvModel = ConductivityTransform->PhysicalToGeneralized(InvModel);

    std::ofstream startmodfile("start.mod");
    std::copy(InvModel.begin(), InvModel.end(), std::ostream_iterator<double>(
        startmodfile, "\n"));
    boost::shared_ptr<jiba::X3DObjective>
        X3DObjective(new jiba::X3DObjective());
    X3DObjective->SetObservedData(Data);
    X3DObjective->SetCoarseModelGeometry(Model);
    X3DObjective->SetDataCovar(MTDataError);

    boost::shared_ptr<jiba::JointObjective> Objective(new jiba::JointObjective(
        true));
    boost::shared_ptr<jiba::MatOpRegularization> Regularization =
        RegSetup.SetupObjective(vm, Model, ConductivityTransform, CovModVec);

    double lambda = 1.0;
    std::cout << "Lambda: ";
    std::cin >> lambda;
    Objective->AddObjective(X3DObjective, ConductivityTransform, 1.0, "MT");
    Objective->AddObjective(Regularization, ConductivityTransform, lambda,
        "Regularization");

    size_t maxiter = 30;
    std::cout << "Maximum number of iterations: ";
    std::cin >> maxiter;
    std::cout << "Performing inversion." << std::endl;

    boost::shared_ptr<jiba::GradientBasedOptimization> Optimizer =
        InversionSetup.ConfigureInversion(vm, Objective, InvModel, CovModVec);

    std::ofstream misfitfile("misfit.out");
    misfitfile << "0 " << Objective->CalcMisfit(InvModel) << " ";
    std::ofstream difffile("diff.out");
    std::copy(X3DObjective->GetDataDifference().begin(),
        X3DObjective->GetDataDifference().end(), std::ostream_iterator<double>(
            difffile, "\n"));
    std::flush(difffile);

    std::copy(Objective->GetIndividualFits().begin(),
        Objective->GetIndividualFits().end(), std::ostream_iterator<double>(
            misfitfile, " "));
    misfitfile << std::endl;

    size_t iteration = 0;
    boost::posix_time::ptime starttime =
        boost::posix_time::microsec_clock::local_time();

    bool terminate = false;
    do
      {
        std::cout << "Iteration: " << iteration + 1 << std::endl;
        Optimizer->MakeStep(InvModel);
        std::copy(InvModel.begin(), InvModel.end(),
            Model.SetConductivities().origin());
        Model.WriteVTK(modelfilename + jiba::stringify(iteration)
            + ".mt.raw.vtk");

        jiba::rvec CondInvModel = ConductivityTransform->GeneralizedToPhysical(
            InvModel);
        std::copy(CondInvModel.begin(), CondInvModel.end(),
            Model.SetConductivities().origin());
        Model.WriteVTK(modelfilename + jiba::stringify(iteration)
            + ".mt.inv.vtk");
        std::cout << "Currrent Misfit: " << Optimizer->GetMisfit() << std::endl;
        std::cout << "Currrent Gradient: " << Optimizer->GetGradNorm()
            << std::endl;
        std::cout << std::endl;

        misfitfile << iteration + 1 << " " << Optimizer->GetMisfit() << " ";
        std::copy(Objective->GetIndividualFits().begin(),
            Objective->GetIndividualFits().end(),
            std::ostream_iterator<double>(misfitfile, " "));
        misfitfile << " " << Objective->GetNEval();
        misfitfile << std::endl;
        iteration++;
        terminate = jiba::WantAbort();
      } while (iteration < maxiter && !terminate
        && Objective->GetIndividualFits()[0] > ndata);

    boost::posix_time::ptime endtime =
        boost::posix_time::microsec_clock::local_time();
    double cachedruntime = (endtime - starttime).total_seconds();
    std::cout << "Runtime: " << cachedruntime << " s" << std::endl;
    std::cout << std::endl;
    InvModel = ConductivityTransform->GeneralizedToPhysical(InvModel);

    std::copy(InvModel.begin(), InvModel.begin()
        + Model.GetConductivities().num_elements(),
        Model.SetConductivities().origin());

    //calculate the predicted data
    std::cout << "Calculating response of inversion model." << std::endl;
    jiba::rvec InvData(jiba::X3DMTCalculator().Calculate(Model));
    jiba::WriteImpedancesToNetCDF(modelfilename + ".inv_imp.nc", Frequencies,
        XCoord, YCoord, ZCoord, InvData);

    //and write out the data and model
    //here we have to distinguish again between scalar and ftg data
    std::cout << "Writing out inversion results." << std::endl;

    Model.WriteVTK(modelfilename + ".inv.vtk");
    Model.WriteNetCDF(modelfilename + ".inv.nc");
    std::cout << std::endl;
  }
