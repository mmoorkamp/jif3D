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
#ifdef HAVEOPENMP
#include <omp.h>
#endif
#include <boost/program_options.hpp>
#include "../Global/convert.h"
#include "../Global/FatalException.h"
#include "../Global/NumUtil.h"
#include "../Global/VectorTransform.h"
#include "../Global/FileUtil.h"
#include "../Global/Noise.h"
#include "../ModelBase/VTKTools.h"
#include "../ModelBase/NetCDFModelTools.h"
#include "../Inversion/LimitedMemoryQuasiNewton.h"
#include "../Inversion/NonLinearConjugateGradient.h"
#include "../Inversion/JointObjective.h"
#include "../Regularization/MinDiffRegularization.h"
#include "../Inversion/ModelTransforms.h"
#include "../Gravity/ReadWriteGravityData.h"
#include "../GravMag/ThreeDGravMagCalculator.h"
#include "../GravMag/MinMemGravMagCalculator.h"
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
    desc.add_options()("help", "produce help message")("threads", po::value<int>(),
        "The number of openmp threads")("covmod", po::value<std::string>(),
        "A file containing the model covariance");

    jif3D::SetupRegularization RegSetup;
    jif3D::SetupInversion InversionSetup;
    jif3D::SetupGravity GravitySetup;
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
#ifdef HAVEOPENMP
    if (vm.count("threads"))
      {
        omp_set_num_threads(vm["threads"].as<int>());
      }
#endif
    jif3D::rvec CovModVec;
    if (vm.count("covmod"))
      {
        jif3D::ThreeDGravityModel CovModel;
        CovModel.ReadNetCDF(vm["covmod"].as<std::string>());
        const size_t ncovmod = CovModel.GetDensities().num_elements();
        CovModVec.resize(ncovmod);
        std::copy(CovModel.GetDensities().origin(),
            CovModel.GetDensities().origin() + ncovmod, CovModVec.begin());
      }

    boost::shared_ptr<jif3D::GeneralModelTransform> Transform;

    boost::shared_ptr<jif3D::JointObjective> Objective(new jif3D::JointObjective(true));
    GravitySetup.SetupObjective(vm, *Objective.get(), Transform);

    jif3D::rvec InvModel(GravitySetup.GetScalModel().GetDensities().num_elements());
    std::copy(GravitySetup.GetScalModel().GetDensities().origin(),
        GravitySetup.GetScalModel().GetDensities().origin()
            + GravitySetup.GetScalModel().GetDensities().num_elements(),
        InvModel.begin());

    boost::shared_ptr<jif3D::ObjectiveFunction> Regularization = RegSetup.SetupObjective(
        vm, GravitySetup.GetScalModel(), CovModVec);

    double reglambda = 1.0;

    std::cout << "Regularization Lambda: ";
    std::cin >> reglambda;
    Objective->AddObjective(Regularization, Transform, reglambda);

    std::cout << "Performing inversion." << std::endl;

    //jif3D::rvec Ones(CovModVec);
    //std::fill(Ones.begin(),Ones.end(),1.0);
    boost::shared_ptr<jif3D::GradientBasedOptimization> Optimizer =
        InversionSetup.ConfigureInversion(vm, Objective, InvModel, CovModVec);

    size_t iteration = 0;
    size_t maxiter = 30;
    std::cout << "Maximum number of iterations: ";
    std::cin >> maxiter;
    std::string modelfilename("result");
    std::ofstream misfitfile("misfit.out");
    jif3D::ThreeDGravityModel GravModel(GravitySetup.GetScalModel());
    do
      {
        try
          {
            std::cout << "Iteration: " << iteration << std::endl;
            Optimizer->MakeStep(InvModel);

            jif3D::rvec DensInvModel = Transform->GeneralizedToPhysical(InvModel);
            std::copy(DensInvModel.begin(), DensInvModel.end(),
                GravModel.SetDensities().origin());
            GravModel.WriteVTK(
                modelfilename + jif3D::stringify(iteration) + ".grav.inv.vtk");
            GravModel.WriteNetCDF(
                modelfilename + jif3D::stringify(iteration) + ".grav.inv.nc");

            ++iteration;
            std::cout << "Gradient Norm: " << Optimizer->GetGradNorm() << std::endl;
            std::cout << "Currrent Misfit: " << Optimizer->GetMisfit() << std::endl;
            std::cout << "Currrent Gradient: " << Optimizer->GetGradNorm() << std::endl;
            misfitfile << std::setw(5) << iteration << " " << std::setw(15)
                << Optimizer->GetMisfit() << " ";
            for (size_t i = 0; i < Objective->GetIndividualFits().size(); ++i)
              {
                misfitfile << std::setw(15) << Objective->GetIndividualFits().at(i)
                    << " ";
              }

            misfitfile << " " << Objective->GetNEval();
            misfitfile << std::endl;
          } catch (jif3D::FatalException &e)
          {
            std::cerr << e.what() << std::endl;
            iteration = maxiter;
          }
      } while (iteration < maxiter && Optimizer->GetMisfit() > 1
        && Optimizer->GetGradNorm() > 1e-6);

    jif3D::rvec DensInvModel(Transform->GeneralizedToPhysical(InvModel));

    std::copy(DensInvModel.begin(),
        DensInvModel.begin() + GravModel.SetDensities().num_elements(),
        GravModel.SetDensities().origin());

    //calculate the predicted data
    std::cout << "Calculating response of inversion model." << std::endl;

    typedef typename jif3D::MinMemGravMagCalculator<jif3D::ThreeDGravityModel> CalculatorType;
    boost::shared_ptr<CalculatorType> ScalGravityCalculator = boost::shared_ptr<
        CalculatorType>(jif3D::CreateGravityCalculator<CalculatorType>::MakeScalar());
    jif3D::rvec GravInvData(ScalGravityCalculator->Calculate(GravModel));
    jif3D::SaveScalarGravityMeasurements(modelfilename + ".inv_sgd.nc", GravInvData,
        GravModel.GetMeasPosX(), GravModel.GetMeasPosY(), GravModel.GetMeasPosZ(),
        GravitySetup.GetScalGravObjective().GetDataError());

    boost::shared_ptr<CalculatorType> FTGGravityCalculator = boost::shared_ptr<
        CalculatorType>(jif3D::CreateGravityCalculator<CalculatorType>::MakeTensor());
    jif3D::rvec FTGInvData(FTGGravityCalculator->Calculate(GravModel));
    jif3D::SaveTensorGravityMeasurements(modelfilename + ".inv_ftg.nc", FTGInvData,
        GravModel.GetMeasPosX(), GravModel.GetMeasPosY(), GravModel.GetMeasPosZ(),
        GravitySetup.GetFTGObjective().GetDataError());

    jif3D::Write3DDataToVTK(modelfilename + ".inv_sgd.vtk", "grav_accel", GravInvData,
        GravModel.GetMeasPosX(), GravModel.GetMeasPosY(), GravModel.GetMeasPosZ());
    jif3D::Write3DTensorDataToVTK(modelfilename + ".inv_ftg.vtk", "U", FTGInvData,
        GravModel.GetMeasPosX(), GravModel.GetMeasPosY(), GravModel.GetMeasPosZ());
    //and write out the data and model
    //here we have to distinguish again between scalar and ftg data
    std::cout << "Writing out inversion results." << std::endl;
    GravModel.WriteVTK(modelfilename + ".grav.inv.vtk");
    GravModel.WriteNetCDF(modelfilename + ".grav.inv.nc");
    std::cout << std::endl;
  }
