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
#include <boost/mpi.hpp>
#include <boost/bind.hpp>
#include <boost/program_options.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>

#include <boost/mpi/packed_oarchive.hpp>
#include <boost/mpi/packed_iarchive.hpp>
#include <boost/serialization/shared_ptr.hpp>
#include <boost/serialization/serialization.hpp>
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
#include "../Inversion/MPIThreeDModelObjective.h"
#include "../MT/X3DModel.h"
#include "../MT/X3DMTCalculator.h"
#include "../MT/ReadWriteImpedances.h"
#include "../Joint/SetupRegularization.h"
#include "../Joint/SetupInversion.h"

namespace ublas = boost::numeric::ublas;
namespace po = boost::program_options;


//BOOST_CLASS_EXPORT( jif3D::VectorTransform )
BOOST_CLASS_EXPORT( jif3D::CopyTransform )

BOOST_CLASS_EXPORT( jif3D::GeneralModelTransform )
BOOST_CLASS_EXPORT( jif3D::ModelCopyTransform )
BOOST_CLASS_EXPORT( jif3D::TanhTransform )
BOOST_CLASS_EXPORT( jif3D::DensityTransform )
//BOOST_CLASS_EXPORT( jif3D::ChainedTransform )
//BOOST_CLASS_EXPORT( jif3D::SectionTransform )
//BOOST_CLASS_EXPORT( jif3D::MultiSectionTransform )
BOOST_CLASS_EXPORT( jif3D::ObjectiveFunction )
BOOST_CLASS_EXPORT( jif3D::X3DMTCalculator )
BOOST_CLASS_EXPORT( jif3D::X3DModel )
BOOST_CLASS_EXPORT( jif3D::ModelRefiner )
//BOOST_CLASS_EXPORT( jif3D::ThreeDModelObjective<jif3D::X3DMTCalculator> )
//BOOST_CLASS_EXPORT( jif3D::MPIThreeDModelObjective<jif3D::X3DMTCalculator> )
//BOOST_CLASS_EXPORT( jif3D::JointObjective )
int main(int argc, char *argv[])
  {


    double mincond = 1e-6;
    double maxcond = 10;

    boost::mpi::environment env(argc, argv);
    boost::mpi::communicator world;
    bool is_generator = world.rank() > 0;
    boost::mpi::communicator local = world.split(is_generator ? 0 : 1);
    std::vector<int> mpiindices(local.size());
    for (int i = 0; i < world.size(); ++i)
      {
        if (i > 0)
          mpiindices.push_back(i);
      }
    boost::shared_ptr<jif3D::JointObjective> Objective(new jif3D::JointObjective(true));

    po::options_description desc("General options");
    desc.add_options()("help", "produce help message")("threads", po::value<int>(),
        "The number of openmp threads")("covmod", po::value<std::string>(),
        "A file containing the model covariance")("mincond",
        po::value(&mincond)->default_value(1e-6),
        "The minimum value for conductivity in S/m")("maxcond",
        po::value(&maxcond)->default_value(10),
        "The maximum value for conductivity in S/m");

    jif3D::SetupRegularization RegSetup;
    jif3D::SetupInversion InversionSetup;
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
        omp_set_num_threads(vm["threads"].as<int>());
      }

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

    //first we read in the starting model and the measured data
    std::string modelfilename = jif3D::AskFilename("Starting model Filename: ");

    //we read in the starting modelfile
    jif3D::X3DModel Model;
    Model.ReadNetCDF(modelfilename);

    //get the name of the file containing the data and read it in
    std::string datafilename = jif3D::AskFilename("Data Filename: ");

    //read in data
    jif3D::rvec Data, ZError;
    std::vector<double> XCoord, YCoord, ZCoord, Frequencies;
    jif3D::ReadImpedancesFromNetCDF(datafilename, Frequencies, XCoord, YCoord, ZCoord,
        Data, ZError);
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
    jif3D::Write3DDataToVTK(datafilename + ".vtk", "Z", FirstFreq, XCoord, YCoord, ZCoord);

    //we define a few constants that are used throughout the inversion
    const size_t ndata = Data.size();
    const double errorlevel = 0.02;

    //create objects for the misfit and a very basic error estimate
    jif3D::rvec MTDataError = jif3D::ConstructMTError(Data, errorlevel);
    std::transform(ZError.begin(), ZError.end(), MTDataError.begin(), ZError.begin(),
        std::max<double>);

    for (size_t i = 0; i < Model.GetConductivities().shape()[2]; ++i)
      {
        Model.SetConductivities()[0][0][i] *= (1 + 0.0001 * (i + 1));
      }
    jif3D::rvec InvModel(Model.GetConductivities().num_elements());
    std::copy(Model.GetConductivities().origin(),
        Model.GetConductivities().origin() + Model.GetConductivities().num_elements(),
        InvModel.begin());
    Model.WriteVTK("start.vtk");

    boost::shared_ptr<jif3D::ChainedTransform> ConductivityTransform(
        new jif3D::ChainedTransform);

    jif3D::rvec RefModel(InvModel);
    std::fill(RefModel.begin(), RefModel.end(), 1.0);
    //because the tanh transform is used inside a logarithmic transform
    //we need to take the natural logarithm of the actual minimum and maximum
    ConductivityTransform->AddTransform(
        boost::shared_ptr<jif3D::GeneralModelTransform>(
            new jif3D::TanhTransform(std::log(mincond), std::log(maxcond))));
    ConductivityTransform->AddTransform(
        boost::shared_ptr<jif3D::GeneralModelTransform>(new jif3D::LogTransform(RefModel)));

    InvModel = ConductivityTransform->PhysicalToGeneralized(InvModel);

    std::ofstream startmodfile("start.mod");
    std::copy(InvModel.begin(), InvModel.end(),
        std::ostream_iterator<double>(startmodfile, "\n"));

    jif3D::X3DMTCalculator Calculator;
    boost::shared_ptr<jif3D::MPIThreeDModelObjective<jif3D::X3DMTCalculator> > X3DObjective(
        new jif3D::MPIThreeDModelObjective<jif3D::X3DMTCalculator>(Calculator, local, world,
            0, mpiindices));

    X3DObjective->SetObservedData(Data);
    X3DObjective->SetCoarseModelGeometry(Model);
    X3DObjective->SetDataError(ZError);

    boost::shared_ptr<jif3D::MatOpRegularization> Regularization = RegSetup.SetupObjective(
        vm, Model, ConductivityTransform, CovModVec);

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

    boost::mpi::broadcast(world, Calculator, 0);

    boost::shared_ptr<jif3D::GradientBasedOptimization> Optimizer =
        InversionSetup.ConfigureInversion(vm, Objective, InvModel, CovModVec);

    std::ofstream misfitfile("misfit.out");
    misfitfile << "0 " << Objective->CalcMisfit(InvModel) << " ";
    std::ofstream difffile("diff.out");
    std::copy(X3DObjective->GetDataDifference().begin(),
        X3DObjective->GetDataDifference().end(),
        std::ostream_iterator<double>(difffile, "\n"));
    std::flush(difffile);

    std::copy(Objective->GetIndividualFits().begin(),
        Objective->GetIndividualFits().end(),
        std::ostream_iterator<double>(misfitfile, " "));
    misfitfile << std::endl;

    size_t iteration = 0;
    boost::posix_time::ptime starttime = boost::posix_time::microsec_clock::local_time();

    bool terminate = false;
    do
      {
        try
          {
            std::cout << "Iteration: " << iteration + 1 << std::endl;
            Optimizer->MakeStep(InvModel);
            std::copy(InvModel.begin(), InvModel.end(),
                Model.SetConductivities().origin());
            Model.WriteVTK(modelfilename + jif3D::stringify(iteration) + ".mt.raw.vtk");

            jif3D::rvec CondInvModel = ConductivityTransform->GeneralizedToPhysical(
                InvModel);
            std::copy(CondInvModel.begin(), CondInvModel.end(),
                Model.SetConductivities().origin());
            Model.WriteVTK(modelfilename + jif3D::stringify(iteration) + ".mt.inv.vtk");
            std::cout << "Current Misfit: " << Optimizer->GetMisfit() << std::endl;
            std::cout << "Current Gradient: " << Optimizer->GetGradNorm() << std::endl;
            std::cout << std::endl;

            misfitfile << iteration + 1 << " " << Optimizer->GetMisfit() << " ";
            std::copy(Objective->GetIndividualFits().begin(),
                Objective->GetIndividualFits().end(),
                std::ostream_iterator<double>(misfitfile, " "));
            misfitfile << " " << Objective->GetNEval();
            misfitfile << std::endl;
            iteration++;
            terminate = jif3D::WantAbort();
          } catch (jif3D::FatalException &e)
          {
            std::cerr << e.what() << std::endl;
            terminate = true;
          }
      } while (iteration < maxiter && !terminate
        && Objective->GetIndividualFits()[0] > ndata);

    boost::posix_time::ptime endtime = boost::posix_time::microsec_clock::local_time();
    double cachedruntime = (endtime - starttime).total_seconds();
    std::cout << "Runtime: " << cachedruntime << " s" << std::endl;
    std::cout << std::endl;
    InvModel = ConductivityTransform->GeneralizedToPhysical(InvModel);

    std::copy(InvModel.begin(),
        InvModel.begin() + Model.GetConductivities().num_elements(),
        Model.SetConductivities().origin());

    //calculate the predicted data
    std::cout << "Calculating response of inversion model." << std::endl;
    jif3D::rvec InvData(jif3D::X3DMTCalculator().Calculate(Model));
    jif3D::WriteImpedancesToNetCDF(modelfilename + ".inv_imp.nc", Frequencies, XCoord,
        YCoord, ZCoord, InvData);
    jif3D::WriteImpedancesToNetCDF(modelfilename + ".diff_imp.nc", Frequencies, XCoord,
        YCoord, ZCoord, X3DObjective->GetDataDifference());

    //and write out the data and model
    //here we have to distinguish again between scalar and ftg data
    std::cout << "Writing out inversion results." << std::endl;

    Model.WriteVTK(modelfilename + ".inv.vtk");
    Model.WriteNetCDF(modelfilename + ".inv.nc");
    std::cout << std::endl;
  }
