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

/*! \file jointinv.cpp
 * The main joint inversion program.
 */

int main(int argc, char *argv[])
  {
    jiba::SetupTomo TomoSetup;
    jiba::SetupGravity GravitySetup;
    jiba::SetupMT MTSetup;
    jiba::SetupInversion InversionSetup;
    jiba::SetupRegularization RegSetup;
    jiba::SetupCoupling CouplingSetup;

    po::options_description desc("General options");
    desc.add_options()("help", "produce help message")("threads",
        po::value<int>(), "The number of openmp threads")("covmod", po::value<
        std::string>(), "A file containing the model covariance");

    desc.add(TomoSetup.SetupOptions());
    desc.add(GravitySetup.SetupOptions());
    desc.add(MTSetup.SetupOptions());
    desc.add(InversionSetup.SetupOptions());
    desc.add(RegSetup.SetupOptions());
    desc.add(CouplingSetup.SetupOptions());

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    //we can also read options from a configuration file
    const std::string ConfFileName("jointinv.conf");
    if (boost::filesystem::exists(ConfFileName))
      {

        std::ifstream ConfFile(ConfFileName.c_str());
        po::store(po::parse_config_file(ConfFile, desc), vm);
      }
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
        jiba::ThreeDSeismicModel CovModel;
        CovModel.ReadNetCDF(vm["covmod"].as<std::string> ());
        const size_t ncovmod = CovModel.GetSlownesses().num_elements();
        CovModVec.resize(ncovmod);
        std::copy(CovModel.GetSlownesses().origin(),
            CovModel.GetSlownesses().origin() + ncovmod, CovModVec.begin());
      }
    //we want some output so we set Verbose in the constructor to true
    boost::shared_ptr<jiba::JointObjective> Objective(new jiba::JointObjective(
        true));

    boost::shared_ptr<jiba::GeneralModelTransform> TomoTransform, MTTransform,
        GravityTransform, RegTransform;

    CouplingSetup.SetupTransforms(vm, TomoTransform, GravityTransform,
        MTTransform, RegTransform);

    jiba::ThreeDSeismicModel TomoModel;
    TomoSetup.SetupObjective(vm, *Objective.get(), TomoModel, TomoTransform);
    GravitySetup.SetupObjective(vm, *Objective.get(), GravityTransform);

    if (!EqualGridGeometry(TomoModel, GravitySetup.GetModel()))
      {
        throw jiba::FatalException(
            "Gravity model does not have the same geometry as starting model");
      }

    MTSetup.SetupObjective(vm, *Objective.get(), MTTransform);
    if (!EqualGridGeometry(MTSetup.GetModel(), TomoModel))
      {
        throw jiba::FatalException(
            "MT model does not have the same geometry as starting model");
      }

    boost::shared_ptr<jiba::MatOpRegularization> Regularization =
        RegSetup.SetupObjective(vm, TomoModel, RegTransform, CovModVec);

    jiba::rvec InvModel;
    CouplingSetup.SetupModelVector(vm, InvModel, TomoModel,
        GravitySetup.GetModel(), MTSetup.GetModel(), *Objective.get(),
        Regularization, RegSetup.GetSubStart());

    size_t maxiter = 1;
    std::cout << "Maximum iterations: ";
    std::cin >> maxiter;

    boost::posix_time::ptime starttime =
        boost::posix_time::microsec_clock::local_time();

    std::cout << "Performing inversion." << std::endl;

    boost::shared_ptr<jiba::GradientBasedOptimization> Optimizer =
        InversionSetup.ConfigureInversion(vm, Objective, InvModel, CovModVec);

    size_t iteration = 0;

    jiba::rvec TomoInvModel(TomoTransform->GeneralizedToPhysical(InvModel));
    std::ofstream misfitfile("misfit.out");
    //calculate initial misfit
    misfitfile << "0 " << Objective->CalcMisfit(InvModel) << " ";
    std::copy(Objective->GetIndividualFits().begin(),
        Objective->GetIndividualFits().end(), std::ostream_iterator<double>(
            misfitfile, " "));
    misfitfile << std::endl;

    std::string modelfilename = "result";
    bool terminate = true;
    do
      {
        terminate = true;
        try
          {
            std::cout << "\n\n Iteration: " << iteration << std::endl;
            //std::copy(InvModel.begin(), InvModel.begin() + ngrid,
            //    TomoModel.SetSlownesses().origin());
            //TomoModel.WriteVTK("raw_model" + jiba::stringify(iteration)
            //    + ".tomo.inv.vtk");
            Optimizer->MakeStep(InvModel);

            ++iteration;

            TomoInvModel = TomoTransform->GeneralizedToPhysical(InvModel);
            std::copy(TomoInvModel.begin(), TomoInvModel.begin()
                + TomoModel.GetNModelElements(),
                TomoModel.SetSlownesses().origin());
            TomoModel.WriteVTK(modelfilename + jiba::stringify(iteration)
                + ".tomo.inv.vtk");
            std::cout << "Currrent Misfit: " << Optimizer->GetMisfit()
                << std::endl;
            std::cout << "Currrent Gradient: " << Optimizer->GetGradNorm()
                << std::endl;
            misfitfile << std::setw(5) << iteration << " " << std::setw(15)
                << Optimizer->GetMisfit() << " ";
            for (size_t i = 0; i < Objective->GetIndividualFits().size(); ++i)
              {
                misfitfile << std::setw(15)
                    << Objective->GetIndividualFits().at(i) << " ";
              }

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
                if (Objective->GetObjective(i).ConvergenceLimit() > 0.0)
                  {
                    std::cout << "Reached target misfit." << std::endl;
                    std::cout << "Objective number: " << i << std::endl;
                    std::cout << "Misfit: "
                        << Objective->GetIndividualFits().at(i) << std::endl;
                    std::cout << "Target: "
                        << Objective->GetObjective(i).GetNData() << std::endl;
                  }
              }
          }
        terminate = jiba::WantAbort();
      } while (iteration < maxiter && !terminate && Optimizer->GetGradNorm()
        > 1e-6);

    jiba::rvec DensInvModel(GravityTransform->GeneralizedToPhysical(InvModel));
    jiba::rvec CondInvModel(MTTransform->GeneralizedToPhysical(InvModel));
    std::copy(TomoInvModel.begin(), TomoInvModel.begin()
        + TomoModel.GetNModelElements(), TomoModel.SetSlownesses().origin());

    jiba::ThreeDGravityModel GravModel(GravitySetup.GetModel());
    std::copy(DensInvModel.begin(), DensInvModel.begin()
        + GravModel.SetDensities().num_elements(),
        GravModel.SetDensities().origin());

    jiba::X3DModel MTModel(MTSetup.GetModel());
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
/* @} */
