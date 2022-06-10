//============================================================================
// Name        : jointinvtomodcmag.cpp
// Author      : August 3, 2014
// Version     :
// Copyright   : 2014, zhanjie
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
#include "../Inversion/ModelTransforms.h"
#include "../Regularization/CrossGradient.h"
#include "../Regularization/MinDiffRegularization.h"
#include "../Tomo/ThreeDSeismicModel.h"
#include "../Tomo/ReadWriteTomographyData.h"
#include "../Tomo/TomographyCalculator.h"
#include "../DCResistivity/ThreeDDCResistivityModel.h"
#include "../DCResistivity/ReadWriteDCResistivityData.h"
#include "../DCResistivity/DCResistivityCalculator.h"
#include "../Magnetics/ThreeDSusceptibilityModel.h"
#include "../Magnetics/ReadWriteMagneticData.h"
#include "../Magnetics/OMPMagneticSusceptibilityImp.h"
#include "../GravMag/ThreeDGravMagImplementation.h"
#include "../GravMag/ThreeDGravMagCalculator.h"
#include "SetupTomo.h"
#include "SetupDCResistivity.h"
#include "SetupMagnetics.h"
#include "SetupInversion.h"
#include "SetupRegularization.h"
#include "SetupTomoDCMagCoupling.h"
#include "InversionOutput.h"

namespace ublas = boost::numeric::ublas;
namespace po = boost::program_options;

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
    jif3D::SetupDCResistivity DCResistivitySetup;
    jif3D::SetupMagnetics MagneticSetup;
    jif3D::SetupInversion InversionSetup;
    jif3D::SetupRegularization RegSetup;
    jif3D::SetupTomoDCMagCoupling CouplingSetup;
    bool WaveletParm = false;
    double xorigin, yorigin;
    double coolingfactor = 1.0;
    //we also create a number of options that are specific to our joint inversion
    //or act globally so that they cannot be associated with one subsystem
    po::options_description desc("General options");
    desc.add_options()("help", "produce help message")("threads", po::value<int>(),
        "The number of openmp threads")("covmod", po::value<std::string>(),
        "A file containing the model covariance")("tempdir", po::value<std::string>(),
        "The name of the directory to store temporary files in")("wavelet",
        "Parametrize inversion by wavelet coefficients")("xorigin",
        po::value(&xorigin)->default_value(0.0),
        "The origin for the inversion grid in x-direction")("yorigin",
        po::value(&yorigin)->default_value(0.0),
        "The origin for the inversion grid in y-direction")("coolingfactor",
        po::value(&coolingfactor)->default_value(1.0),
        "The factor to multiply the weight for the regularization at each iteration EXPERIMENTAL");
    //we need to add the description for each part to the boost program options object
    //that way the user can get a help output and the parser object recongnizes these options
    desc.add(TomoSetup.SetupOptions());
    desc.add(DCResistivitySetup.SetupOptions());
    desc.add(MagneticSetup.SetupOptions());
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
        std::string version = "$Id: jointinv.cpp 600 2014-02-11 10:15:08Z mmoorkamp $";
        std::cout << version << std::endl;
        std::cout << desc << "\n";
        return 1;
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
    boost::shared_ptr<jif3D::GeneralModelTransform> TomoTransform, DCTransform,
        MagTransform;
    //coupling setup is responsible to set the appropriate transformation
    //as we have to use different ones depending on the chose coupling mechanism
    jif3D::ThreeDSeismicModel StartModel;
    CouplingSetup.SetupTransforms(vm, StartModel, TomoTransform, DCTransform,
    		MagTransform, WaveletParm);

    //read in the tomography model and setup the options that are applicable
    //for the seismic tomography part of the inversion
    bool havetomo = TomoSetup.SetupObjective(vm, *Objective.get(), TomoTransform, xorigin,
        yorigin);
    if (havetomo && !EqualGridGeometry(TomoSetup.GetModel(), StartModel))
      {
        throw jif3D::FatalException(
            "Tomography model does not have the same geometry as starting model");
      }
    //setup the dcresistivity part of the joint inversion
    bool havedc = DCResistivitySetup.SetupObjective(vm, *Objective.get(), DCTransform,
        xorigin, yorigin);
    //if we have a seismic and a gravity objective function, we have
    //to make sure that the starting models have the same geometry (not considering refinement)
    if (havedc && !EqualGridGeometry(StartModel, DCResistivitySetup.GetModel()))
      {
        throw jif3D::FatalException(
            "DCResistivity model does not have the same geometry as starting model");
      }
    //setup the MT part of the joint inversion
    bool havemag = MagneticSetup.SetupObjective(vm, *Objective.get(), MagTransform, xorigin,
        yorigin, TempDir);
    //if we have a seismic and a MT objective function, we have
    //to make sure that the starting models have the same geometry (not considering refinement)
    if (havemag && !EqualGridGeometry(MagneticSetup.GetModel(), StartModel))
      {
        throw jif3D::FatalException(
            "Magnetic model does not have the same geometry as starting model");
      }
    //now we setup the regularization
    boost::shared_ptr<jif3D::RegularizationFunction> Regularization =
        RegSetup.SetupObjective(vm, StartModel, CovModVec);

    //the vector InvModel will hold the current inversion model
    //depending on the chosen coupling mechanism it will have different size
    //so we fill its content in the object Coupling setup
    jif3D::rvec InvModel;
    CouplingSetup.SetupModelVector(vm, InvModel, StartModel, TomoSetup.GetModel(),
    		DCResistivitySetup.GetModel(), MagneticSetup.GetModel(), *Objective.get(), Regularization,
        RegSetup.GetSubStart());
    //finally ask for the maximum number of iterations
    size_t maxiter = 1;
    std::cout << "Maximum iterations: ";
    std::cin >> maxiter;
    //note the start time of the core calculations for statistics
    //and output some status information
    boost::posix_time::ptime starttime = boost::posix_time::microsec_clock::local_time();

    std::cout << "Performing inversion." << std::endl;

    auto CovObj = boost::make_shared<jif3D::DiagonalCovariance>(CovModVec);

    boost::shared_ptr<jif3D::GradientBasedOptimization> Optimizer =
        InversionSetup.ConfigureInversion(vm, Objective, InvModel, CovObj);

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


    size_t iteration = 0;
    std::ofstream misfitfile("misfit.out");
    std::ofstream rmsfile("rms.out");
    std::ofstream weightfile("weights.out");
    //calculate initial misfit
    double InitialMisfit = Objective->CalcMisfit(InvModel);
    StoreMisfit(misfitfile, 0, InitialMisfit, *Objective);
    StoreRMS(rmsfile, 0, *Objective);
    StoreWeights(weightfile, 0, *Objective);
    jif3D::ThreeDDCResistivityModel DCModel(DCResistivitySetup.GetModel());
    jif3D::ThreeDSusceptibilityModel MagModel(MagneticSetup.GetModel());
    jif3D::ThreeDSeismicModel TomoModel(TomoSetup.GetModel());

    if (!havetomo)
      TomoModel = StartModel;
    if (!havedc)
      DCModel = StartModel;
    if (!havemag)
      MagModel = StartModel;

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
            SaveModel(InvModel, *DCTransform.get(), DCModel,
                modelfilename + jif3D::stringify(iteration) + ".dc.inv");
            SaveModel(InvModel, *MagTransform.get(), MagModel,
                modelfilename + jif3D::stringify(iteration) + ".mag.inv");
            //write out some information about misfit to the screen
            std::cout << "Currrent Misfit: " << Optimizer->GetMisfit() << std::endl;
            std::cout << "Currrent Gradient: " << Optimizer->GetGradNorm() << std::endl;
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

    SaveModel(InvModel, *TomoTransform.get(), TomoModel, modelfilename + ".tomo.inv");
    SaveModel(InvModel, *DCTransform.get(), DCModel, modelfilename + ".dc.inv");
    SaveModel(InvModel, *MagTransform.get(), MagModel, modelfilename + ".mag.inv");

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
