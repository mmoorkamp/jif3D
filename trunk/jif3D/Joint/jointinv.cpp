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
#ifdef HAVEOPENMP
#include <omp.h>
#endif
#include <boost/program_options.hpp>
#include <boost/program_options/config.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/log/core.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/expressions.hpp>

#include "../Global/convert.h"
#include "../Global/FatalException.h"
#include "../Global/NumUtil.h"
#include "../Global/VectorTransform.h"
#include "../Global/FileUtil.h"
#include "../ModelBase/VTKTools.h"
#include "../ModelBase/NetCDFModelTools.h"
#include "../ModelBase/EqualGeometry.h"
#include "../Inversion/JointObjective.h"
#include "../Inversion/ThreeDModelObjective.h"
#include "../Regularization/MinDiffRegularization.h"
#include "../Inversion/ModelTransforms.h"
#include "../Tomo/ThreeDSeismicModel.h"
#include "../Tomo/ReadWriteTomographyData.h"
#include "../Tomo/TomographyCalculator.h"
#include "../Gravity/ReadWriteGravityData.h"
#include "../Gravity/ThreeDGravityFactory.h"
#include "../MT/X3DModel.h"
#include "../MT/X3DMTCalculator.h"
#include "../MT/ReadWriteImpedances.h"
#include "SetupTomo.h"
#include "SetupGravity.h"
#include "SetupMT.h"
#include "SetupInversion.h"
#include "SetupRegularization.h"
#include "SetupCoupling.h"
#include "InversionOutput.h"

namespace ublas = boost::numeric::ublas;
namespace po = boost::program_options;
namespace logging = boost::log;
/** \addtogroup joint Joint inversion routines */
/* @{ */

/*! \file jointinv.cpp
 * The main joint inversion program. The main task of the program is to read in the appropriate files
 * and options and from these settings assemble the objects that perform the actual work. Also,
 * the main program manages the output of inversion results and data among other statistics.
 */

int main(int argc, char *argv[])
  {
    //first we create objects that manage the setup of individual parts
    //that way we can reuse these objects so that other programs use
    //exactly the same options
    jif3D::SetupTomo TomoSetup;
    jif3D::SetupGravity GravitySetup;
    jif3D::SetupMT MTSetup;
    jif3D::SetupInversion InversionSetup;
    jif3D::SetupRegularization RegSetup;
    jif3D::SetupCoupling CouplingSetup;
    bool WaveletParm = false;
    bool WantSequential = false;
    double xorigin, yorigin;
    double coolingfactor = 1.0;
    //we also create a number of options that are specific to our joint inversion
    //or act globally so that they cannot be associated with one subsystem
    po::options_description desc("General options");
    desc.add_options()("help", "produce help message")("debug",
        "Write debugging information")("threads", po::value<int>(),
        "The number of openmp threads")("covmod", po::value<std::string>(),
        "A file containing the model covariance")("tempdir", po::value<std::string>(),
        "The name of the directory to store temporary files in")("wavelet",
        "Parametrize inversion by wavelet coefficients")("xorigin",
        po::value(&xorigin)->default_value(0.0),
        "The origin for the inversion grid in x-direction")("yorigin",
        po::value(&yorigin)->default_value(0.0),
        "The origin for the inversion grid in y-direction")("coolingfactor",
        po::value(&coolingfactor)->default_value(1.0),
        "The factor to multiply the weight for the regularization at each iteration EXPERIMENTAL")(
        "sequential",
        "Do not create a single objective function, but split into on OF per method");
    //we need to add the description for each part to the boost program options object
    //that way the user can get a help output and the parser object recongnizes these options
    desc.add(TomoSetup.SetupOptions());
    desc.add(GravitySetup.SetupOptions());
    desc.add(MTSetup.SetupOptions());
    desc.add(InversionSetup.SetupOptions());
    desc.add(RegSetup.SetupOptions());
    desc.add(CouplingSetup.SetupOptions());

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    //we can also read options from a configuration file
    //that way we do not have to specify a large number of options
    //on the command line, the order we use now, reading the configuration
    //file after parsing the command line options means that
    //the command line options override configuration file options
    //as one would expect (see also program_options documentation)
    const std::string ConfFileName("jointinv.conf");
    if (boost::filesystem::exists(ConfFileName))
      {
        std::ifstream ConfFile(ConfFileName.c_str());
        po::store(po::parse_config_file(ConfFile, desc), vm);
      }
    po::notify(vm);

    //if the option was "help" we output the program version
    //and a description of all options, but we do not perform any inversion
    if (vm.count("help"))
      {
        std::string version = "$Id$";
        std::cout << version << std::endl;
        std::cout << desc << "\n";
        return 1;
      }

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

    if (vm.count("sequential"))
      {
        WantSequential = true;
      }
    if (vm.count("wavelet"))
      {
        WaveletParm = true;
      }
#ifdef HAVEOPENMP
    if (vm.count("threads"))
      {
        omp_set_num_threads(vm["threads"].as<int>());
      }
#endif
    //some objects accept a directory name as a path to store temporary files
    //we check that this directory actually exists to avoid double checking
    //and early detection of problems
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

    jif3D::rvec CovModVec;
    if (vm.count("covmod"))
      {
        jif3D::ThreeDSeismicModel CovModel;
        //we store the covariances in a seismic model file
        //but we do not have to have an equidistant grid
        CovModel.ReadNetCDF(vm["covmod"].as<std::string>(), false);
        const size_t ncovmod = CovModel.GetSlownesses().num_elements();
        CovModVec.resize(ncovmod);
        std::copy(CovModel.GetSlownesses().origin(),
            CovModel.GetSlownesses().origin() + ncovmod, CovModVec.begin());
      }
    //we want some output so we set Verbose in the constructor to true
    boost::shared_ptr<jif3D::JointObjective> Objective;
    Objective = boost::shared_ptr<jif3D::JointObjective>(new jif3D::JointObjective(true));

    //we need a number of transformation objects to translate the generalized model parameters
    //used by the inversion algorithm to physically meaningful parameters for the forward
    //calculation. The exact type of transformation depends on the chosen coupling
    //and is assigned below
    boost::shared_ptr<jif3D::GeneralModelTransform> TomoTransform, MTTransform,
        GravityTransform;
    //coupling setup is responsible to set the appropriate transformation
    //as we have to use different ones depending on the chose coupling mechanism
    jif3D::ThreeDSeismicModel StartModel;
    CouplingSetup.SetupTransforms(vm, StartModel, TomoTransform, GravityTransform,
        MTTransform, WaveletParm);

    //read in the tomography model and setup the options that are applicable
    //for the seismic tomography part of the inversion
    bool havetomo = TomoSetup.SetupObjective(vm, *Objective.get(), TomoTransform, xorigin,
        yorigin);
    if (havetomo && !EqualGridGeometry(TomoSetup.GetModel(), StartModel))
      {
        throw jif3D::FatalException(
            "Tomography model does not have the same geometry as starting model");
      }
    //setup the gravity part of the joint inversion
    bool havegrav = GravitySetup.SetupObjective(vm, *Objective.get(), GravityTransform,
        xorigin, yorigin, TempDir);
    //if we have a seismic and a gravity objective function, we have
    //to make sure that the starting models have the same geometry (not considering refinement)
    if (havegrav && !EqualGridGeometry(StartModel, GravitySetup.GetScalModel()))
      {
        throw jif3D::FatalException(
            "Gravity model does not have the same geometry as starting model");
      }
    //setup the MT part of the joint inversion
    bool havemt = MTSetup.SetupObjective(vm, *Objective.get(), MTTransform, xorigin,
        yorigin, TempDir);
    //if we have a seismic and a MT objective function, we have
    //to make sure that the starting models have the same geometry (not considering refinement)
    if (havemt && !EqualGridGeometry(MTSetup.GetModel(), StartModel))
      {
        throw jif3D::FatalException(
            "MT model does not have the same geometry as starting model");
      }
    //now we setup the regularization
    boost::shared_ptr<jif3D::RegularizationFunction> Regularization =
        RegSetup.SetupObjective(vm, StartModel, CovModVec);

    //the vector InvModel will hold the current inversion model
    //depending on the chosen coupling mechanism it will have different size
    //so we fill its content in the object Coupling setup
    jif3D::rvec InvModel;
    CouplingSetup.SetupModelVector(vm, InvModel, StartModel, TomoSetup.GetModel(),
        GravitySetup.GetScalModel(), MTSetup.GetModel(), *Objective.get(), Regularization,
        RegSetup.GetSubStart());
    //finally ask for the maximum number of iterations
    size_t maxiter = 1;
    std::cout << "Maximum iterations: ";
    std::cin >> maxiter;
    //note the start time of the core calculations for statistics
    //and output some status information
    boost::posix_time::ptime starttime = boost::posix_time::microsec_clock::local_time();

    std::string modelfilename = "result";
    //write out the seismic source and receiver positions for plotting
    //and general quality control
    if (havetomo)
      {
        jif3D::Write3DDataToVTK(modelfilename + ".rec.vtk", "Receiver",
            jif3D::rvec(TomoSetup.GetModel().GetMeasPosX().size()),
            TomoSetup.GetModel().GetMeasPosX(), TomoSetup.GetModel().GetMeasPosY(),
            TomoSetup.GetModel().GetMeasPosZ());
        jif3D::Write3DDataToVTK(modelfilename + ".sor.vtk", "Source",
            jif3D::rvec(TomoSetup.GetModel().GetSourcePosX().size()),
            TomoSetup.GetModel().GetSourcePosX(), TomoSetup.GetModel().GetSourcePosY(),
            TomoSetup.GetModel().GetSourcePosZ());
      }

    if (havemt)
      {
        jif3D::Write3DDataToVTK(modelfilename + ".mt_sites.vtk", "MT Sites",
            jif3D::rvec(MTSetup.GetModel().GetMeasPosX().size()),
            MTSetup.GetModel().GetMeasPosX(), MTSetup.GetModel().GetMeasPosY(),
            MTSetup.GetModel().GetMeasPosZ());
      }

    std::cout << "Calculating initial misfit." << std::endl;
    size_t iteration = 0;
    std::ofstream misfitfile("misfit.out");
    std::ofstream rmsfile("rms.out");
    std::ofstream weightfile("weights.out");
    //calculate initial misfit
    double InitialMisfit = Objective->CalcMisfit(InvModel);
    StoreMisfit(misfitfile, 0, InitialMisfit, *Objective);
    StoreRMS(rmsfile, 0, *Objective);
    StoreWeights(weightfile, 0, *Objective);
    jif3D::ThreeDGravityModel GravModel(GravitySetup.GetScalModel());
    jif3D::X3DModel MTModel(MTSetup.GetModel());
    jif3D::ThreeDSeismicModel TomoModel(TomoSetup.GetModel());

    if (!havetomo)
      TomoModel = StartModel;
    if (!havemt)
      MTModel = StartModel;
    if (!havegrav)
      GravModel = StartModel;

    std::cout << "Performing inversion." << std::endl;
    if (WantSequential)
      {
        const size_t ngrid = StartModel.GetSlownesses().num_elements();
        boost::shared_ptr<jif3D::JointObjective> TomoObjective(Objective->clone());
        boost::shared_ptr<jif3D::JointObjective> MTObjective(Objective->clone());
        boost::shared_ptr<jif3D::JointObjective> GravObjective(Objective->clone());

        std::ofstream tomomisfitfile("tomo_misfit.out");
        std::ofstream mtmisfitfile("mt_misfit.out");
        std::ofstream gravmisfitfile("grav_misfit.out");
        std::vector<double> Weights = Objective->GetWeights();
        std::vector<double> TomoWeights(Weights);
        std::vector<double> MTWeights(Weights);
        std::vector<double> GravWeights(Weights);
        TomoWeights.at(1) = 0.0;
        TomoWeights.at(2) = 0.0;
        TomoWeights.at(3) = 0.0;
        TomoWeights.at(6) = 0.0;
        TomoWeights.at(8) = 0.0;
        TomoWeights.at(9) = 0.0;
        MTWeights.at(0) = 0.0;
        MTWeights.at(1) = 0.0;
        MTWeights.at(2) = 0.0;
        MTWeights.at(4) = 0.0;
        MTWeights.at(7) = 0.0;
        MTWeights.at(8) = 0.0;
        GravWeights.at(0) = 0.0;
        GravWeights.at(3) = 0.0;
        GravWeights.at(5) = 0.0;
        GravWeights.at(7) = 0.0;
        GravWeights.at(9) = 0.0;

        TomoObjective->SetWeights(TomoWeights);
        MTObjective->SetWeights(MTWeights);
        GravObjective->SetWeights(GravWeights);

        jif3D::rvec TomoInvModel = InvModel;
        jif3D::rvec MTInvModel = InvModel;
        jif3D::rvec GravInvModel = InvModel;

        boost::shared_ptr<jif3D::GradientBasedOptimization> TomoOptimizer =
            InversionSetup.ConfigureInversion(vm, TomoObjective, TomoInvModel, CovModVec);
        boost::shared_ptr<jif3D::GradientBasedOptimization> MTOptimizer =
            InversionSetup.ConfigureInversion(vm, MTObjective, MTInvModel, CovModVec);
        boost::shared_ptr<jif3D::GradientBasedOptimization> GravOptimizer =
            InversionSetup.ConfigureInversion(vm, GravObjective, GravInvModel, CovModVec);

        bool terminate = false;
        jif3D::rvec OldModel(InvModel);
        //this is the core inversion loop, we make optimization steps
        //until either we reach the maximum number of iterations
        //or fulfill a termination criterion
        while (iteration < maxiter && !terminate)
          {
            terminate = true;
            jif3D::rvec OldTomo = TomoInvModel;
            jif3D::rvec OldGrav = GravInvModel;
            jif3D::rvec OldMT = MTInvModel;

            //we catch all jif3D internal exceptions so that we can graciously
            //exit and write out some final information before stopping the program

            std::cout << "\n\n Iteration: " << iteration << std::endl;

            //update the inversion model
            try
              {
                TomoOptimizer->MakeStep(TomoInvModel);
              } catch (jif3D::FatalException &e)
              {
                TomoInvModel = OldTomo;
                std::cerr << "In tomography: " << e.what() << std::endl;
              }
            try
              {
                MTOptimizer->MakeStep(MTInvModel);
              } catch (jif3D::FatalException &e)
              {
                MTInvModel = OldMT;
                std::cerr << "In MT: " << e.what() << std::endl;
              }
            try
              {
                GravOptimizer->MakeStep(GravInvModel);
              } catch (jif3D::FatalException &e)
              {
                GravInvModel = OldGrav;
                std::cerr << "In gravity: " << e.what() << std::endl;
              }

            ublas::subrange(MTInvModel, 0, ngrid) = ublas::subrange(TomoInvModel, 0,
                ngrid);
            ublas::subrange(GravInvModel, 0, ngrid) = ublas::subrange(TomoInvModel, 0,
                ngrid);

            ublas::subrange(TomoInvModel, ngrid, 2 * ngrid) = ublas::subrange(
                GravInvModel, ngrid, 2 * ngrid);
            ublas::subrange(MTInvModel, ngrid, 2 * ngrid) = ublas::subrange(GravInvModel,
                ngrid, 2 * ngrid);

            ublas::subrange(TomoInvModel, 2 * ngrid, 3 * ngrid) = ublas::subrange(
                MTInvModel, 2 * ngrid, 3 * ngrid);
            ublas::subrange(GravInvModel, 2 * ngrid, 3 * ngrid) = ublas::subrange(
                MTInvModel, 2 * ngrid, 3 * ngrid);
            //Objective->MultiplyWeights(jif3D::JointObjective::regularization,
            //    coolingfactor);
            ++iteration;
            //we save all models at each iteration, so we can look at the development
            // and use intermediate models in case something goes wrong
            SaveModel(TomoInvModel, *TomoTransform.get(), TomoModel,
                modelfilename + jif3D::stringify(iteration) + ".tomo.inv");
            SaveModel(MTInvModel, *MTTransform.get(), MTModel,
                modelfilename + jif3D::stringify(iteration) + ".mt.inv");
            SaveModel(GravInvModel, *GravityTransform.get(), GravModel,
                modelfilename + jif3D::stringify(iteration) + ".grav.inv");

            //and write the current misfit for all objectives to a misfit file
            StoreMisfit(tomomisfitfile, iteration, TomoOptimizer->GetMisfit(),
                *TomoObjective);
            StoreMisfit(mtmisfitfile, iteration, MTOptimizer->GetMisfit(), *MTObjective);
            StoreMisfit(gravmisfitfile, iteration, GravOptimizer->GetMisfit(),
                *GravObjective);

            //StoreRMS(rmsfile, iteration, *Objective);
            //StoreWeights(weightfile, iteration, *Objective);
            std::cout << "\n\n";
            //we stop when either we do not make any improvement any more
            terminate = CheckConvergence(*Objective);
            //or the file abort exists in the current directory
            terminate = terminate || jif3D::WantAbort();
          }

        ublas::subrange(InvModel, 0, ngrid) = ublas::subrange(TomoInvModel, 0, ngrid);
        ublas::subrange(InvModel, ngrid, 2 * ngrid) = ublas::subrange(GravInvModel, ngrid,
            2 * ngrid);
        ublas::subrange(InvModel, 2 * ngrid, 3 * ngrid) = ublas::subrange(MTInvModel,
            2 * ngrid, 3 * ngrid);
      }
    else
      {

        boost::shared_ptr<jif3D::GradientBasedOptimization> Optimizer =
            InversionSetup.ConfigureInversion(vm, Objective, InvModel, CovModVec);

        bool terminate = false;
        jif3D::rvec OldModel(InvModel);
        //this is the core inversion loop, we make optimization steps
        //until either we reach the maximum number of iterations
        //or fulfill a termination criterion
        while (iteration < maxiter && !terminate)
          {
            terminate = true;
            //we catch all jif3D internal exceptions so that we can graciously
            //exit and write out some final information before stopping the program
            try
              {
                std::cout << "\n\n Iteration: " << iteration << std::endl;
                //we save the current model so we can go back to it
                //in case the optimization step fails
                OldModel = InvModel;
                //update the inversion model
                Optimizer->MakeStep(InvModel);

                Objective->MultiplyWeights(jif3D::JointObjective::regularization,
                    coolingfactor);
                ++iteration;
                //we save all models at each iteration, so we can look at the development
                // and use intermediate models in case something goes wrong
                SaveModel(InvModel, *TomoTransform.get(), TomoModel,
                    modelfilename + jif3D::stringify(iteration) + ".tomo.inv");
                SaveModel(InvModel, *MTTransform.get(), MTModel,
                    modelfilename + jif3D::stringify(iteration) + ".mt.inv");
                SaveModel(InvModel, *GravityTransform.get(), GravModel,
                    modelfilename + jif3D::stringify(iteration) + ".grav.inv");
                //write out some information about misfit to the screen
                std::cout << "Currrent Misfit: " << Optimizer->GetMisfit() << std::endl;
                std::cout << "Currrent Gradient: " << Optimizer->GetGradNorm()
                    << std::endl;
                //and write the current misfit for all objectives to a misfit file
                StoreMisfit(misfitfile, iteration, Optimizer->GetMisfit(), *Objective);
                StoreRMS(rmsfile, iteration, *Objective);
                StoreWeights(weightfile, iteration, *Objective);
                std::cout << "\n\n";
              } catch (jif3D::FatalException &e)
              {
                std::cerr << e.what() << std::endl;
                InvModel = OldModel;
                iteration = maxiter;
              }
            //we stop when either we do not make any improvement any more
            terminate = CheckConvergence(*Objective);
            //or the file abort exists in the current directory
            terminate = terminate || jif3D::WantAbort();
          }
      }
    SaveModel(InvModel, *TomoTransform.get(), TomoModel, modelfilename + ".tomo.inv");
    SaveModel(InvModel, *MTTransform.get(), MTModel, modelfilename + ".mt.inv");
    SaveModel(InvModel, *GravityTransform.get(), GravModel, modelfilename + ".grav.inv");

    //calculate the predicted refraction data
    std::cout << "Calculating response of inversion model." << std::endl;
    //during the last iteration we might have performed steps in the line search
    //so we update the forward calculation for the last proper inversion model
    //we use the results from this calculation to save the final inversion synthetic data
    if (iteration > 0)
      {
        Objective->CalcMisfit(InvModel);
      }
    if (havetomo)
      {
        jif3D::rvec TomoInvData(TomoSetup.GetTomoObjective().GetSyntheticData());
        jif3D::rvec TomoError(TomoSetup.GetTomoObjective().GetDataError());
        jif3D::SaveTraveltimes(modelfilename + ".inv_tt.nc", TomoInvData, TomoError,
            TomoModel);
        jif3D::rvec TomoDiff(TomoSetup.GetTomoObjective().GetIndividualMisfit());
        jif3D::SaveTraveltimes(modelfilename + ".diff_tt.nc", TomoDiff, TomoError,
            TomoModel);
      }
    //if we are inverting gravity data and have specified site locations
    if (havegrav)
      {
        //and write out the data and model
        //here we have to distinguish again between scalar and ftg data
        if (GravitySetup.GetHaveScal())
          {
            jif3D::rvec ScalGravInvData(
                GravitySetup.GetScalGravObjective().GetSyntheticData());
            jif3D::SaveScalarGravityMeasurements(modelfilename + ".inv_sgd.nc",
                ScalGravInvData, GravModel.GetMeasPosX(), GravModel.GetMeasPosY(),
                GravModel.GetMeasPosZ(),
                GravitySetup.GetScalGravObjective().GetDataError());
            jif3D::rvec ScalDiff(
                GravitySetup.GetScalGravObjective().GetIndividualMisfit());
            jif3D::SaveScalarGravityMeasurements(modelfilename + ".diff_sgd.nc", ScalDiff,
                GravModel.GetMeasPosX(), GravModel.GetMeasPosY(), GravModel.GetMeasPosZ(),
                GravitySetup.GetScalGravObjective().GetDataError());
            jif3D::Write3DDataToVTK(modelfilename + ".inv_sgd.vtk", "grav_accel",
                ScalGravInvData, GravModel.GetMeasPosX(), GravModel.GetMeasPosY(),
                GravModel.GetMeasPosZ());

          }
        if (GravitySetup.GetHaveFTG())
          {
            jif3D::rvec FTGInvData(GravitySetup.GetFTGObjective().GetSyntheticData());
            jif3D::SaveTensorGravityMeasurements(modelfilename + ".inv_ftg.nc",
                FTGInvData, GravModel.GetMeasPosX(), GravModel.GetMeasPosY(),
                GravModel.GetMeasPosZ(), GravitySetup.GetFTGObjective().GetDataError());
            jif3D::rvec FTGDiff(GravitySetup.GetFTGObjective().GetIndividualMisfit());
            jif3D::SaveTensorGravityMeasurements(modelfilename + ".diff_ftg.nc", FTGDiff,
                GravModel.GetMeasPosX(), GravModel.GetMeasPosY(), GravModel.GetMeasPosZ(),
                GravitySetup.GetFTGObjective().GetDataError());
            jif3D::Write3DTensorDataToVTK(modelfilename + ".diff_ftg.vtk", "U",
                FTGInvData, GravModel.GetMeasPosX(), GravModel.GetMeasPosY(),
                GravModel.GetMeasPosZ());
          }
      }

    //if we are inverting MT data and have specified site locations
    if (havemt)
      {
        std::cout << "C: ";
        std::copy(MTModel.GetDistortionParameters().begin(),
            MTModel.GetDistortionParameters().end(),
            std::ostream_iterator<double>(std::cout, "\n"));
        //calculate MT inversion result
        jif3D::rvec MTInvData(MTSetup.GetMTObjective().GetSyntheticData());
        jif3D::WriteImpedancesToNetCDF(modelfilename + ".inv_mt.nc",
            MTModel.GetFrequencies(), MTModel.GetMeasPosX(), MTModel.GetMeasPosY(),
            MTModel.GetMeasPosZ(), MTInvData, MTSetup.GetMTObjective().GetDataError(),
            MTModel.GetDistortionParameters());
        jif3D::rvec MTDiff(MTSetup.GetMTObjective().GetIndividualMisfit());
        jif3D::WriteImpedancesToNetCDF(modelfilename + ".diff_mt.nc",
            MTModel.GetFrequencies(), MTModel.GetMeasPosX(), MTModel.GetMeasPosY(),
            MTModel.GetMeasPosZ(), MTDiff, MTSetup.GetMTObjective().GetDataError());
      }

    std::ofstream datadiffile("data.diff");
    std::copy(Objective->GetDataDifference().begin(),
        Objective->GetDataDifference().end(),
        std::ostream_iterator<double>(datadiffile, "\n"));
    boost::posix_time::ptime endtime = boost::posix_time::microsec_clock::local_time();
    double cachedruntime = (endtime - starttime).total_seconds();
    std::cout << "Runtime: " << cachedruntime << " s" << std::endl;
    std::cout << std::endl;
  }
/* @} */
