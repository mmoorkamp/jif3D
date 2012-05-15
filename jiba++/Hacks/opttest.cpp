//============================================================================
// Name        : opttest.cpp
// Author      : Sep 19, 2008
// Version     :
// Copyright   : 2008, mmoorkamp
//============================================================================

#ifdef HAVE_CONFIG_H
#include "OPT++_config.h"
#endif

#include <fstream>
#include <limits>
#ifdef HAVE_STD
#include <cstdio>
#else
#include <stdio.h>
#endif
#include <numeric>
#include <cmath>
#include <omp.h>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/program_options.hpp>
#include <boost/bind.hpp>
#include "OptppDecl.h"
#include "../Global/convert.h"
#include "../Global/FatalException.h"
#include "../Global/NumUtil.h"
#include "../Global/VectorTransform.h"
#include "../Global/FileUtil.h"
#include "../Global/Noise.h"
#include "../Inversion/ModelTransforms.h"
#include "../Inversion/JointObjective.h"
#include "../Regularization/MinDiffRegularization.h"
#include "../Regularization/GradientRegularization.h"
#include "../ModelBase/VTKTools.h"
#include "../ModelBase/NetCDFTools.h"
#include "../ModelBase/EqualGeometry.h"
#include "../Gravity/ReadWriteGravityData.h"
#include "../Gravity/ThreeDGravityCalculator.h"
#include "../Gravity/MinMemGravityCalculator.h"
#include "../Gravity/ThreeDGravityFactory.h"
#include "../Tomo/ThreeDSeismicModel.h"
#include "../Tomo/ReadWriteTomographyData.h"
#include "../Tomo/TomographyCalculator.h"
#include "../Joint/SetupTomo.h"
#include "../Joint/SetupGravity.h"
#include "../Joint/SetupMT.h"
#include "../Joint/SetupInversion.h"
#include "../Joint/SetupRegularization.h"
#include "../Joint/SetupCoupling.h"

using NEWMAT::ColumnVector;
namespace po = boost::program_options;

jiba::JointObjective &GetObjective()
  {
    // see  http://www.parashift.com/c%2B%2B-faq-lite/ctors.html#faq-10.15
    //why we have to do it like this.
    static jiba::JointObjective *Objective = new jiba::JointObjective(true);
    return *Objective;
  }

jiba::rvec InvModel;

void update_model(int, int, ColumnVector x)
  {
    for (size_t i = 0; i < InvModel.size(); ++i)
      InvModel(i) = x(i + 1);
  }

void init_model(int ndim, ColumnVector& x)
  {
    if (ndim != InvModel.size())
      {
        std::cerr << "Ndim: " << ndim << " does not match Model size " << InvModel.size()
            << std::endl;
        exit(1);
      }
    for (size_t i = 0; i < InvModel.size(); ++i)
      x(i + 1) = InvModel(i);
  }

void eval_objective(int mode, int n, const ColumnVector& x, double& fx, ColumnVector& g,
    int& result)
  {

    if (n != InvModel.size())
      {
        std::cerr << "N: " << n << " does not match Model size " << InvModel.size()
            << std::endl;
        return;
      }
    for (size_t i = 0; i < InvModel.size(); ++i)
      InvModel(i) = x(i + 1);

    if (mode & OPTPP::NLPFunction)
      {
        fx = GetObjective().CalcMisfit(InvModel);
        result = OPTPP::NLPFunction;
      }
    if (mode & OPTPP::NLPGradient)
      {
        GetObjective().CalcMisfit(InvModel);
        jiba::rvec Gradient(GetObjective().CalcGradient(InvModel));
        for (size_t i = 0; i < Gradient.size(); ++i)
          g(i + 1) = Gradient(i);

        result = OPTPP::NLPGradient;
      }
  }

int main(int argc, char *argv[])
  {

    jiba::SetupTomo TomoSetup;
    jiba::SetupGravity GravitySetup;
    jiba::SetupMT MTSetup;
    jiba::SetupRegularization RegSetup;
    jiba::SetupCoupling CouplingSetup;

    po::options_description desc("General options");
    desc.add_options()("help", "produce help message")("threads", po::value<int>(),
        "The number of openmp threads")("covmod", po::value<std::string>(),
        "A file containing the model covariance");

    desc.add(TomoSetup.SetupOptions());
    desc.add(GravitySetup.SetupOptions());
    desc.add(MTSetup.SetupOptions());
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
        omp_set_num_threads(vm["threads"].as<int>());
      }

    jiba::rvec CovModVec;
    if (vm.count("covmod"))
      {
        jiba::ThreeDSeismicModel CovModel;
        //we store the covariances in a seismic model file
        //but we do not have to have an equidistant grid
        CovModel.ReadNetCDF(vm["covmod"].as<std::string>(), false);
        const size_t ncovmod = CovModel.GetSlownesses().num_elements();
        CovModVec.resize(ncovmod);
        std::copy(CovModel.GetSlownesses().origin(),
            CovModel.GetSlownesses().origin() + ncovmod, CovModVec.begin());
      }

    boost::shared_ptr<jiba::GeneralModelTransform> TomoTransform, MTTransform,
        GravityTransform, RegTransform;

    jiba::ThreeDSeismicModel StartModel;
    CouplingSetup.SetupTransforms(vm, StartModel, TomoTransform, GravityTransform,
        MTTransform, RegTransform);

    bool havetomo = TomoSetup.SetupObjective(vm, GetObjective(), TomoTransform);
    bool havegrav = GravitySetup.SetupObjective(vm, GetObjective(), GravityTransform);

    if (havetomo && havegrav && !EqualGridGeometry(StartModel, GravitySetup.GetScalModel()))
      {
        throw jiba::FatalException(
            "Gravity model does not have the same geometry as starting model");
      }

    bool havemt = MTSetup.SetupObjective(vm, GetObjective(), MTTransform);
    if (havetomo && havemt && !EqualGridGeometry(MTSetup.GetModel(), StartModel))
      {
        throw jiba::FatalException(
            "MT model does not have the same geometry as starting model");
      }

    boost::shared_ptr<jiba::MatOpRegularization> Regularization = RegSetup.SetupObjective(
        vm, StartModel, RegTransform, CovModVec);
    CouplingSetup.SetupModelVector(vm, InvModel, StartModel, GravitySetup.GetScalModel(),
        MTSetup.GetModel(), GetObjective(), Regularization, RegSetup.GetSubStart());

    size_t maxiter = 1;
    std::cout << "Maximum iterations: ";
    std::cin >> maxiter;

    boost::posix_time::ptime starttime = boost::posix_time::microsec_clock::local_time();

    std::cout << "Performing inversion." << std::endl;

    size_t iteration = 0;
    std::ofstream misfitfile("misfit.out");
    //calculate initial misfit
    misfitfile << "0 " << GetObjective().CalcMisfit(InvModel) << " ";
    std::copy(GetObjective().GetIndividualFits().begin(),
        GetObjective().GetIndividualFits().end(),
        std::ostream_iterator<double>(misfitfile, " "));
    misfitfile << std::endl;

    std::string modelfilename = "result";

    jiba::Write3DDataToVTK(modelfilename + ".rec.vtk", "Receiver",
        jiba::rvec(TomoSetup.GetModel().GetMeasPosX().size()),
        TomoSetup.GetModel().GetMeasPosX(), TomoSetup.GetModel().GetMeasPosY(),
        TomoSetup.GetModel().GetMeasPosZ());
    jiba::Write3DDataToVTK(modelfilename + ".sor.vtk", "Source",
        jiba::rvec(TomoSetup.GetModel().GetSourcePosX().size()),
        TomoSetup.GetModel().GetSourcePosX(), TomoSetup.GetModel().GetSourcePosY(),
        TomoSetup.GetModel().GetSourcePosZ());

    jiba::ThreeDGravityModel GravModel(GravitySetup.GetScalModel());
    jiba::X3DModel MTModel(MTSetup.GetModel());

    if (!havemt)
      MTModel = StartModel;
    if (!havegrav)
      GravModel = StartModel;

    std::cout << "Performing inversion." << std::endl;

    static char *status_file =
      { "tstLBFGS.out" };
    int n = InvModel.size();
    //  Create a Nonlinear problem object

    OPTPP::NLF1 nlp(n, eval_objective, init_model);
    nlp.setIsExpensive(1);

    //  Build a LBFGS object and optimize

    OPTPP::OptLBFGS objfcn(&nlp);
    objfcn.setMaxIter(maxiter);
    objfcn.setUpdateModel(update_model);
    if (!objfcn.setOutputFile(status_file, 0))
      cerr << "main: output file open failed" << endl;
    objfcn.setGradTol(1.e-6);
    objfcn.setMaxBacktrackIter(10);
    objfcn.setPrintFinalX(true);
    objfcn.optimize();

    objfcn.printStatus("Solution from LBFGS: More and Thuente's linesearch");

    objfcn.cleanup();

    jiba::rvec DensInvModel(GravityTransform->GeneralizedToPhysical(InvModel));
    jiba::rvec TomoInvModel(TomoTransform->GeneralizedToPhysical(InvModel));

    std::copy(TomoInvModel.begin(),
        TomoInvModel.begin() + StartModel.GetSlownesses().num_elements(),
        StartModel.SetSlownesses().origin());
    std::copy(DensInvModel.begin(),
        DensInvModel.begin() + GravModel.SetDensities().num_elements(),
        GravModel.SetDensities().origin());

    //calculate the predicted refraction data
    std::cout << "Calculating response of inversion model." << std::endl;
    jiba::rvec TomoInvData(jiba::TomographyCalculator().Calculate(StartModel));
    jiba::SaveTraveltimes(modelfilename + ".inv_tt.nc", TomoInvData, StartModel);

    jiba::rvec ScalGravInvData(
        jiba::CreateGravityCalculator<jiba::MinMemGravityCalculator>::MakeScalar()->Calculate(
            GravModel));
    jiba::rvec FTGInvData(
        jiba::CreateGravityCalculator<jiba::MinMemGravityCalculator>::MakeTensor()->Calculate(
            GravModel));
    jiba::SaveScalarGravityMeasurements(modelfilename + ".inv_sgd.nc", ScalGravInvData,
        GravModel.GetMeasPosX(), GravModel.GetMeasPosY(), GravModel.GetMeasPosZ());
    jiba::SaveTensorGravityMeasurements(modelfilename + ".inv_ftg.nc", FTGInvData,
        GravModel.GetMeasPosX(), GravModel.GetMeasPosY(), GravModel.GetMeasPosZ());
    //and write out the data and model
    //here we have to distinguish again between scalar and ftg data
    std::cout << "Writing out inversion results." << std::endl;

    StartModel.WriteNetCDF(modelfilename + ".tomo.inv.nc");
    StartModel.WriteVTK(modelfilename + ".tomo.inv.vtk");
    GravModel.WriteVTK(modelfilename + ".grav.inv.vtk");
    GravModel.WriteNetCDF(modelfilename + ".grav.inv.nc");
    std::cout << "Number of evaluations: " << GetObjective().GetNEval() << std::endl;
    boost::posix_time::ptime endtime = boost::posix_time::microsec_clock::local_time();
    double cachedruntime = (endtime - starttime).total_seconds();
    std::cout << "Runtime: " << cachedruntime << " s" << std::endl;
    std::cout << std::endl;

  }
