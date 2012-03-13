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
#include <boost/program_options.hpp>
#include <boost/program_options/config.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include "../Global/convert.h"
#include "../Global/FatalException.h"
#include "../Global/NumUtil.h"
#include "../Global/VectorTransform.h"
#include "../Global/FileUtil.h"
#include "../ModelBase/VTKTools.h"
#include "../ModelBase/NetCDFTools.h"
#include "../ModelBase/EqualGeometry.h"
#include "../Inversion/JointObjective.h"
#include "../Inversion/ThreeDModelObjective.h"
#include "../Regularization/MinDiffRegularization.h"
#include "../Inversion/ModelTransforms.h"
#include "../Tomo/ThreeDSeismicModel.h"
#include "../Tomo/ReadWriteTomographyData.h"
#include "../Tomo/TomographyCalculator.h"
#include "../Gravity/ReadWriteGravityData.h"
#include "../Gravity/MinMemGravityCalculator.h"
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

namespace ublas = boost::numeric::ublas;
namespace po = boost::program_options;

/** \addtogroup joint Joint inversion routines */
/* @{ */

//! Use a parameter transform to translate the model vector to a model object and write the model to files for plotting and saving
template<class ModelType>
void SaveModel(const jiba::rvec &InvModel,
    const jiba::GeneralModelTransform &Transform, ModelType &ModelObject,
    const std::string &filename)
  {
    jiba::rvec TransModel = Transform.GeneralizedToPhysical(InvModel);
    assert(TransModel.size() >= ModelObject.GetNModelElements());
    std::copy(TransModel.begin(), TransModel.begin()
        + ModelObject.GetNModelElements(), ModelObject.SetData().origin());
    ModelObject.WriteVTK(filename + ".vtk");
    ModelObject.WriteNetCDF(filename + ".nc");
  }

//! Store the current misfit for all individual objectives with appropriate formating in an output stream
void StoreMisfit(std::ofstream &misfitfile, const size_t iteration,
    const double Misfit, const jiba::JointObjective &Objective)
  {
    misfitfile << std::setw(5) << iteration << " " << std::setw(15) << Misfit
        << " ";
    for (size_t i = 0; i < Objective.GetIndividualFits().size(); ++i)
      {
        misfitfile << std::setw(15) << Objective.GetIndividualFits().at(i)
            << " ";
      }

    misfitfile << " " << Objective.GetNEval();
    misfitfile << std::endl;
  }

//! Check whether we have reached the target misfit for one of the objective functions in the JointObjective object
bool CheckConvergence(const jiba::JointObjective &Objective)
  {
    bool terminate = true;
    for (size_t i = 0; i < Objective.GetIndividualFits().size() - 1; ++i)
      {
        if (Objective.GetIndividualFits().at(i)
            > Objective.GetObjective(i).GetNData())
          {
            terminate = false;
          }
        else
          {
            if (Objective.GetObjective(i).ConvergenceLimit() > 0.0)
              {
                std::cout << "Reached target misfit." << std::endl;
                std::cout << "Objective number: " << i << std::endl;
                std::cout << "Misfit: " << Objective.GetIndividualFits().at(i)
                    << std::endl;
                std::cout << "Target: " << Objective.GetObjective(i).GetNData()
                    << std::endl;
              }
          }
      }
    return terminate;
  }

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
    jiba::SetupTomo TomoSetup;
    jiba::SetupGravity GravitySetup;
    jiba::SetupMT MTSetup;
    jiba::SetupInversion InversionSetup;
    jiba::SetupRegularization RegSetup;
    jiba::SetupCoupling CouplingSetup;
    bool WaveletParm = false;
    //we also create a number of options that are specific to our joint inversion
    //or act globally so that they cannot be associated with one subsystem
    po::options_description desc("General options");
    desc.add_options()("help", "produce help message")("threads",
        po::value<int>(), "The number of openmp threads")("covmod", po::value<
        std::string>(), "A file containing the model covariance")("tempdir",
        po::value<std::string>(),
        "The name of the directory to store temporary files in")("wavelet",
        po::value(&WaveletParm)->default_value(false),
        "Parametrize inversion by wavelet coefficients");
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
        std::string version =
            "$Id$";
        std::cout << version << std::endl;
        std::cout << desc << "\n";
        return 1;
      }

    if (vm.count("threads"))
      {
        omp_set_num_threads(vm["threads"].as<int> ());
      }
    //some objects accept a directory name as a path to store temporary files
    //we check that this directory actually exists to avoid double checking
    //and early detection of problems
    boost::filesystem::path TempDir = boost::filesystem::current_path();
    if (vm.count("tempdir"))
      {
        TempDir = vm["tempdir"].as<std::string> ();
        if (!boost::filesystem::is_directory(TempDir))
          {
            std::cerr << TempDir.string()
                << " is not a directory or does not exist ! \n";
            return 500;
          }
      }

    jiba::rvec CovModVec;
    if (vm.count("covmod"))
      {
        jiba::ThreeDSeismicModel CovModel;
        //we store the covariances in a seismic model file
        //but we do not have to have an equidistant grid
        CovModel.ReadNetCDF(vm["covmod"].as<std::string> (), false);
        const size_t ncovmod = CovModel.GetSlownesses().num_elements();
        CovModVec.resize(ncovmod);
        std::copy(CovModel.GetSlownesses().origin(),
            CovModel.GetSlownesses().origin() + ncovmod, CovModVec.begin());
      }
    //we want some output so we set Verbose in the constructor to true
    boost::shared_ptr<jiba::JointObjective> Objective(new jiba::JointObjective(
        true));
    //we need a number of transformation objects to translate the generalized model parameters
    //used by the inversion algorithm to physically meaningful parameters for the forward
    //calculation. The exact type of transformation depends on the chosen coupling
    //and is assigned below
    boost::shared_ptr<jiba::GeneralModelTransform> TomoTransform, MTTransform,
        GravityTransform, RegTransform;
    //coupling setup is responsible to set the appropriate transformation
    //as we have to use different ones depending on the chose coupling mechanism
    CouplingSetup.SetupTransforms(vm, TomoTransform, GravityTransform,
        MTTransform, RegTransform, WaveletParm);

    //read in the tomography model and setup the options that are applicable
    //for the seismic tomography part of the inversion
    jiba::ThreeDSeismicModel TomoModel;
    bool havetomo = TomoSetup.SetupObjective(vm, *Objective.get(), TomoModel,
        TomoTransform);
    //setup the gravity part of the joint inversion
    bool havegrav = GravitySetup.SetupObjective(vm, *Objective.get(),
        GravityTransform, TempDir);
    //if we have a seismic and a gravity objective function, we have
    //to make sure that the starting models have the same geometry (not considering refinement)
    if (havetomo && havegrav && !EqualGridGeometry(TomoModel,
        GravitySetup.GetModel()))
      {
        throw jiba::FatalException(
            "Gravity model does not have the same geometry as starting model");
      }
    //setup the MT part of the joint inversion
    bool havemt = MTSetup.SetupObjective(vm, *Objective.get(), MTTransform,
        TempDir);
    //if we have a seismic and a MT objective function, we have
    //to make sure that the starting models have the same geometry (not considering refinement)
    if (havetomo && havemt && !EqualGridGeometry(MTSetup.GetModel(), TomoModel))
      {
        throw jiba::FatalException(
            "MT model does not have the same geometry as starting model");
      }
    //now we setup the regularization
    boost::shared_ptr<jiba::MatOpRegularization> Regularization =
        RegSetup.SetupObjective(vm, TomoModel, RegTransform, CovModVec);

    //the vector InvModel will hold the current inversion model
    //depending on the chosen coupling mechanism it will have different size
    //so we fill its content in the object Coupling setup
    jiba::rvec InvModel;
    CouplingSetup.SetupModelVector(vm, InvModel, TomoModel,
        GravitySetup.GetModel(), MTSetup.GetModel(), *Objective.get(),
        Regularization, RegSetup.GetSubStart());
    //finally ask for the maximum number of iterations
    size_t maxiter = 1;
    std::cout << "Maximum iterations: ";
    std::cin >> maxiter;
    //note the start time of the core calculations for statistics
    //and output some status information
    boost::posix_time::ptime starttime =
        boost::posix_time::microsec_clock::local_time();

    std::cout << "Performing inversion." << std::endl;

    boost::shared_ptr<jiba::GradientBasedOptimization> Optimizer =
        InversionSetup.ConfigureInversion(vm, Objective, InvModel, CovModVec);

    size_t iteration = 0;
    std::ofstream misfitfile("misfit.out");
    //calculate initial misfit
    double InitialMisfit = Objective->CalcMisfit(InvModel);
    StoreMisfit(misfitfile, 0, InitialMisfit, *Objective);

    std::string modelfilename = "result";
    //write out the seismic source and receiver positions for plotting
    //and general quality control
    jiba::Write3DDataToVTK(modelfilename + ".rec.vtk", "Receiver", jiba::rvec(
        TomoModel.GetMeasPosX().size()), TomoModel.GetMeasPosX(),
        TomoModel.GetMeasPosY(), TomoModel.GetMeasPosZ());
    jiba::Write3DDataToVTK(modelfilename + ".sor.vtk", "Source", jiba::rvec(
        TomoModel.GetSourcePosX().size()), TomoModel.GetSourcePosX(),
        TomoModel.GetSourcePosY(), TomoModel.GetSourcePosZ());

    jiba::ThreeDGravityModel GravModel(GravitySetup.GetModel());
    jiba::X3DModel MTModel(MTSetup.GetModel());

    if (!havemt)
      MTModel = TomoModel;
    if (!havegrav)
      GravModel = TomoModel;

    bool terminate = false;
    jiba::rvec OldModel(InvModel);
    //this is the core inversion loop, we make optimization steps
    //until either we reach the maximum number of iterations
    //or fulfill a termination criterion
    while (iteration < maxiter && !terminate)
      {
        terminate = true;
        //we catch all jiba internal exceptions so that we can graciously
        //exit and write out some final information before stopping the program
        try
          {
            std::cout << "\n\n Iteration: " << iteration << std::endl;
            //we save the current model so we can go back to it
            //in case the optimization step fails
            OldModel = InvModel;
            //update the inversion model
            Optimizer->MakeStep(InvModel);

            ++iteration;
            //we save all models at each iteration, so we can look at the development
            // and use intermediate models in case something goes wrong
            SaveModel(InvModel, *TomoTransform.get(), TomoModel, modelfilename
                + jiba::stringify(iteration) + ".tomo.inv");
            SaveModel(InvModel, *MTTransform.get(), MTModel, modelfilename
                + jiba::stringify(iteration) + ".mt.inv");
            SaveModel(InvModel, *GravityTransform.get(), GravModel,
                modelfilename + jiba::stringify(iteration) + ".grav.inv");
            //write out some information about misfit to the screen
            std::cout << "Currrent Misfit: " << Optimizer->GetMisfit()
                << std::endl;
            std::cout << "Currrent Gradient: " << Optimizer->GetGradNorm()
                << std::endl;
            //and write the current misfit for all objectives to a misfit file
            StoreMisfit(misfitfile, iteration, Optimizer->GetMisfit(),
                *Objective);
            std::cout << "\n\n";
          } catch (jiba::FatalException &e)
          {
            std::cerr << e.what() << std::endl;
            InvModel = OldModel;
            iteration = maxiter;
          }
        //we stop when either we do not make any improvement any more
        terminate = CheckConvergence(*Objective);
        //or the file abort exists in the current directory
        terminate = terminate || jiba::WantAbort();
      }

    SaveModel(InvModel, *TomoTransform.get(), TomoModel, modelfilename
        + ".tomo.inv");
    SaveModel(InvModel, *MTTransform.get(), MTModel, modelfilename + ".mt.inv");
    SaveModel(InvModel, *GravityTransform.get(), GravModel, modelfilename
        + ".grav.inv");

    //calculate the predicted refraction data
    std::cout << "Calculating response of inversion model." << std::endl;
    //during the last iteration we might have performed steps in the line search
    //so we update the forward calculation for the last proper inversion model
    //we use the results from this calculation to save the final inversion synthetic data
    Objective->CalcMisfit(InvModel);

    jiba::rvec TomoInvData(TomoSetup.GetTomoObjective().GetSyntheticData());
    jiba::SaveTraveltimes(modelfilename + ".inv_tt.nc", TomoInvData, TomoModel);

    //if we are inverting gravity data and have specified site locations
    if (havegrav)
      {
        //and write out the data and model
        //here we have to distinguish again between scalar and ftg data
        if (GravitySetup.GetHaveScal())
          {
            jiba::rvec ScalGravInvData(
                GravitySetup.GetScalGravObjective().GetSyntheticData());
            jiba::SaveScalarGravityMeasurements(modelfilename + ".inv_sgd.nc",
                ScalGravInvData, GravModel.GetMeasPosX(),
                GravModel.GetMeasPosY(), GravModel.GetMeasPosZ());
            jiba::Write3DDataToVTK(modelfilename + ".inv_sgd.vtk",
                "grav_accel", ScalGravInvData, GravModel.GetMeasPosX(),
                GravModel.GetMeasPosY(), GravModel.GetMeasPosZ());

          }
        if (GravitySetup.GetHaveFTG())
          {
            jiba::rvec FTGInvData(
                GravitySetup.GetFTGObjective().GetSyntheticData());
            jiba::SaveTensorGravityMeasurements(modelfilename + ".inv_ftg.nc",
                FTGInvData, GravModel.GetMeasPosX(), GravModel.GetMeasPosY(),
                GravModel.GetMeasPosZ());
            jiba::Write3DTensorDataToVTK(modelfilename + ".inv_ftg.vtk", "U",
                FTGInvData, GravModel.GetMeasPosX(), GravModel.GetMeasPosY(),
                GravModel.GetMeasPosZ());
          }
      }

    //if we are inverting MT data and have specified site locations
    if (havemt)
      {
        //calculate MT inversion result
        jiba::rvec MTInvData(MTSetup.GetMTObjective().GetSyntheticData());
        jiba::WriteImpedancesToNetCDF(modelfilename + "inv_mt.nc",
            MTModel.GetFrequencies(), MTModel.GetMeasPosX(),
            MTModel.GetMeasPosY(), MTModel.GetMeasPosZ(), MTInvData,
            MTSetup.GetMTObjective().GetDataCovar());
      }

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
/* @} */
