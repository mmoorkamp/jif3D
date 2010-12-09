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
#include "../Inversion/ModelTransforms.h"
#include "../Gravity/ReadWriteGravityData.h"
#include "../Gravity/ThreeDGravityCalculator.h"
#include "../Gravity/MinMemGravityCalculator.h"
#include "../Gravity/DepthWeighting.h"
#include "../Gravity/ThreeDGravityFactory.h"
#include "../Joint/SetupRegularization.h"
#include "../Joint/SetupInversion.h"
#include "../Joint/SetupGravity.h"

namespace ublas = boost::numeric::ublas;
namespace po = boost::program_options;

int main(int argc, char *argv[])
  {

    po::options_description desc("General options");
    desc.add_options()("help", "produce help message")("threads",
        po::value<int>(), "The number of openmp threads")("covmod", po::value<
        std::string>(), "A file containing the model covariance");

    jiba::SetupRegularization RegSetup;
    jiba::SetupInversion InversionSetup;
    jiba::SetupGravity GravitySetup;
    desc.add(RegSetup.SetupOptions());
    desc.add(InversionSetup.SetupOptions());
    desc.add(GravitySetup.SetupOptions());
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
        jiba::ThreeDGravityModel CovModel;
        CovModel.ReadNetCDF(vm["covmod"].as<std::string> ());
        const size_t ncovmod = CovModel.GetDensities().num_elements();
        CovModVec.resize(ncovmod);
        std::copy(CovModel.GetDensities().origin(),
            CovModel.GetDensities().origin() + ncovmod, CovModVec.begin());
      }

    boost::shared_ptr<jiba::GeneralModelTransform> Transform;

    boost::shared_ptr<jiba::JointObjective> Objective(new jiba::JointObjective(
        true));
    GravitySetup.SetupObjective(vm, *Objective.get(), Transform);

    jiba::rvec InvModel(GravitySetup.GetModel().GetDensities().num_elements());
    std::copy(GravitySetup.GetModel().GetDensities().origin(),
        GravitySetup.GetModel().GetDensities().origin()
            + GravitySetup.GetModel().GetDensities().num_elements(),
        InvModel.begin());



    boost::shared_ptr<jiba::MatOpRegularization> Regularization =
        RegSetup.SetupObjective(vm, GravitySetup.GetModel(), Transform,  CovModVec);

    double scalgravlambda = 1.0;
    double ftglambda = 1.0;
    double reglambda = 1.0;


    std::cout << "Regularization Lambda: ";
    std::cin >> reglambda;
    Objective->AddObjective(Regularization, Transform, reglambda);

    std::cout << "Performing inversion." << std::endl;

    //jiba::rvec Ones(CovModVec);
    //std::fill(Ones.begin(),Ones.end(),1.0);
    boost::shared_ptr<jiba::GradientBasedOptimization> Optimizer =
        InversionSetup.ConfigureInversion(vm, Objective, InvModel, CovModVec);

    size_t iteration = 0;
    size_t maxiter = 30;
    std::cout << "Maximum number of iterations: ";
    std::cin >> maxiter;
    std::string modelfilename("result");
    std::ofstream misfitfile("misfit.out");
    jiba::ThreeDGravityModel GravModel(GravitySetup.GetModel());
    do
      {
        try
          {
            std::cout << "Iteration: " << iteration << std::endl;
            Optimizer->MakeStep(InvModel);

            jiba::rvec DensInvModel =
                Transform->GeneralizedToPhysical(InvModel);
            std::copy(DensInvModel.begin(), DensInvModel.end(),
                GravModel.SetDensities().origin());
            GravModel.WriteVTK(modelfilename + jiba::stringify(iteration)
                + ".grav.inv.vtk");

            ++iteration;
            std::cout << "Gradient Norm: " << Optimizer->GetGradNorm()
                << std::endl;
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
          } catch (jiba::FatalException &e)
          {
            std::cerr << e.what() << std::endl;
            iteration = maxiter;
          }
      } while (iteration < maxiter && Optimizer->GetMisfit() > 1
        && Optimizer->GetGradNorm() > 1e-6);

    jiba::rvec DensInvModel(Transform->GeneralizedToPhysical(InvModel));

    std::copy(DensInvModel.begin(), DensInvModel.begin()
        + GravModel.SetDensities().num_elements(),
        GravModel.SetDensities().origin());

    //calculate the predicted data
    std::cout << "Calculating response of inversion model." << std::endl;

    boost::shared_ptr<jiba::MinMemGravityCalculator>
        ScalGravityCalculator =
            boost::shared_ptr<jiba::MinMemGravityCalculator>(
                jiba::CreateGravityCalculator<jiba::MinMemGravityCalculator>::MakeScalar());
    jiba::rvec GravInvData(ScalGravityCalculator->Calculate(GravModel));
    jiba::SaveScalarGravityMeasurements(modelfilename + ".inv_sgd.nc",
        GravInvData, GravModel.GetMeasPosX(), GravModel.GetMeasPosY(),
        GravModel.GetMeasPosZ());

    boost::shared_ptr<jiba::MinMemGravityCalculator>
        FTGGravityCalculator =
            boost::shared_ptr<jiba::MinMemGravityCalculator>(
                jiba::CreateGravityCalculator<jiba::MinMemGravityCalculator>::MakeTensor());
    jiba::rvec FTGInvData(FTGGravityCalculator->Calculate(GravModel));
    jiba::SaveTensorGravityMeasurements(modelfilename + ".inv_ftg.nc",
        FTGInvData, GravModel.GetMeasPosX(), GravModel.GetMeasPosY(),
        GravModel.GetMeasPosZ());

    jiba::Write3DDataToVTK(modelfilename + ".inv_sgd.vtk", "grav_accel",
        GravInvData, GravModel.GetMeasPosX(), GravModel.GetMeasPosY(),
        GravModel.GetMeasPosZ());
    jiba::Write3DTensorDataToVTK(modelfilename + ".inv_ftg.vtk", "U",
        FTGInvData, GravModel.GetMeasPosX(), GravModel.GetMeasPosY(),
        GravModel.GetMeasPosZ());
    //and write out the data and model
    //here we have to distinguish again between scalar and ftg data
    std::cout << "Writing out inversion results." << std::endl;
    GravModel.WriteVTK(modelfilename + ".grav.inv.vtk");
    GravModel.WriteNetCDF(modelfilename + ".grav.inv.nc");
    std::cout << std::endl;
  }
