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
#include "../Inversion/ModelTransforms.h"
#include "../Inversion/JointObjective.h"
#include "../Inversion/MinDiffRegularization.h"
#include "../Inversion/GradientRegularization.h"
#include "../Inversion/ConstructError.h"
#include "../Inversion/JointObjective.h"
#include "../ModelBase/VTKTools.h"
#include "../ModelBase/NetCDFTools.h"
#include "../Gravity/GravityObjective.h"
#include "../Gravity/ReadWriteGravityData.h"
#include "../Gravity/ThreeDGravityCalculator.h"
#include "../Gravity/MinMemGravityCalculator.h"
#include "../Gravity/DepthWeighting.h"
#include "../Gravity/ThreeDGravityFactory.h"
#include "../Tomo/ThreeDSeismicModel.h"
#include "../Tomo/ReadWriteTomographyData.h"
#include "../Tomo/TomographyObjective.h"
#include "../Tomo/TomographyCalculator.h"


using NEWMAT::ColumnVector;
namespace po = boost::program_options;
boost::shared_ptr<jiba::JointObjective> Objective(new jiba::JointObjective());
jiba::rvec InvModel;

void update_model(int, int, ColumnVector x)
  {
    for (size_t i = 0; i < InvModel.size(); ++i)
      InvModel( i) = x(i + 1);
  }

void init_model(int ndim, ColumnVector& x)
  {
    if (ndim != InvModel.size())
      {
        std::cerr << "Ndim: " << ndim << " does not match Model size "
            << InvModel.size() << std::endl;
        exit(1);
      }
    for (size_t i = 0; i < InvModel.size(); ++i)
      x(i + 1) = InvModel(i);
  }

void eval_objective(int mode, int n, const ColumnVector& x, double& fx,
    ColumnVector& g, int& result)
  {

    if (n != InvModel.size())
      {
        std::cerr << "N: " << n << " does not match Model size "
            << InvModel.size() << std::endl;
        return;
      }
    for (size_t i = 0; i < InvModel.size(); ++i)
      InvModel( i) = x(i + 1);

    if (mode & OPTPP::NLPFunction)
      {
        fx = Objective->CalcMisfit(InvModel);
        result = OPTPP::NLPFunction;
      }
    if (mode & OPTPP::NLPGradient)
      {
        Objective->CalcMisfit(InvModel);
        jiba::rvec Gradient(Objective->CalcGradient());
        for (size_t i = 0; i < Gradient.size(); ++i)
          g(i + 1) = Gradient(i);

        result = OPTPP::NLPGradient;
      }
  }

int main(int argc, char *argv[])
  {

    po::options_description desc("Allowed options");
    desc.add_options()("help", "produce help message")("cpu",
        "Perform calculation on CPU [default]")("gpu",
        "Perform calculation on GPU")("threads", po::value<int>(),
        "The number of openmp threads")("corrpairs", po::value<int>(),
        "The number correction pairs for L-BFGS")("nlcg",
        "Use NLCG optimization");

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);

    if (vm.count("help"))
      {
        std::cout << desc << "\n";
        return 1;
      }
    bool wantcuda = false;

    if (vm.count("gpu"))
      {
        std::cout << "Using GPU" << "\n";
        wantcuda = true;
      }
    if (vm.count("threads"))
      {
        omp_set_num_threads(vm["threads"].as<int> ());
      }
    int correctionpairs = 5;
    if (vm.count("corrpairs"))
      {
        correctionpairs = vm["corrpairs"].as<int> ();
      }
    bool wantnlcg = false;
    if (vm.count("nlcg"))
      {
        wantnlcg = true;
      }
    //these objects hold information about the measurements and their geometry
    jiba::rvec TomoData, ScalGravData, FTGData;

    //first we read in the starting model and the measured data
    std::string modelfilename = jiba::AskFilename("Starting model Filename: ");
    //we read in the starting modelfile
    jiba::ThreeDSeismicModel TomoModel;
    TomoModel.ReadNetCDF(modelfilename);
    TomoModel.WriteVTK(modelfilename + ".vtk");
    //get the name of the file containing the data and read it in
    std::string tomodatafilename = jiba::AskFilename(
        "Tomography Data Filename: ");

    //read in data
    jiba::ReadTraveltimes(tomodatafilename, TomoData, TomoModel);

    std::string scalgravdatafilename = jiba::AskFilename(
        "Scalar Gravity Data Filename: ");
    std::string ftgdatafilename = jiba::AskFilename("FTG Data Filename: ");
    std::string gravmodelfilename = jiba::AskFilename(
        "Gravity Model Filename: ");
    jiba::ThreeDGravityModel GravModel;
    GravModel.ReadNetCDF(gravmodelfilename);
    GravModel = TomoModel;

    jiba::ThreeDGravityModel::tMeasPosVec PosX, PosY, PosZ;
    jiba::ReadScalarGravityMeasurements(scalgravdatafilename, ScalGravData,
        PosX, PosY, PosZ);
    jiba::ReadTensorGravityMeasurements(ftgdatafilename, FTGData, PosX, PosY,
        PosZ);
    GravModel.ClearMeasurementPoints();
    for (size_t i = 0; i < PosX.size(); ++i)
      {
        GravModel.AddMeasurementPoint(PosX.at(i), PosY.at(i), PosZ.at(i));
      }

    //if we don't have data inversion doesn't make sense;
    if (TomoData.empty() || ScalGravData.empty())
      {
        std::cerr << "No measurements defined" << std::endl;
        exit(100);
      }

    InvModel.resize(TomoModel.GetSlownesses().num_elements());
    std::copy(TomoModel.GetSlownesses().origin(),
        TomoModel.GetSlownesses().origin()
            + TomoModel.GetSlownesses().num_elements(), InvModel.begin());

    const double z0 = 5.0;
    const double DepthExponent = -2.0;
    jiba::rvec WeightVector, ModelWeight(InvModel.size());
    //calculate the depth scaling
    jiba::ConstructDepthWeighting(GravModel.GetZCellSizes(), z0, WeightVector,
        jiba::WeightingTerm(DepthExponent));
    for (size_t i = 0; i < ModelWeight.size(); ++i)
      {
        ModelWeight( i) = WeightVector(i % GravModel.GetZCellSizes().size());
      }

    jiba::rvec RefModel(InvModel);
    jiba::rvec PreCond(InvModel.size());
    std::fill(PreCond.begin(), PreCond.end(), 1.0);
    const double minslow = 1e-4;
    const double maxslow = 0.005;
    boost::shared_ptr<jiba::GeneralModelTransform> DensityTransform(
        new jiba::TanhDensityTransform(minslow, maxslow));
    boost::shared_ptr<jiba::GeneralModelTransform> SlownessTransform(
        new jiba::TanhTransform(minslow, maxslow));
    //double average = std::accumulate(InvModel.begin(),InvModel.end(),0.0)/InvModel.size();
    //std::fill(RefModel.begin(),RefModel.end(),average);
    InvModel = SlownessTransform->PhysicalToGeneralized(InvModel);
    jiba::rvec
        DensStartModel(DensityTransform->GeneralizedToPhysical(InvModel));
    std::cout << "Background layers: "
        << GravModel.GetBackgroundDensities().size() << std::endl;
    std::copy(DensStartModel.begin(), DensStartModel.end(),
        GravModel.SetDensities().origin());
    GravModel.WriteNetCDF("out_dens.nc");

    boost::shared_ptr<jiba::TomographyObjective> TomoObjective(
        new jiba::TomographyObjective());
    TomoObjective->SetObservedData(TomoData);
    TomoObjective->SetModelGeometry(TomoModel);
    TomoObjective->SetPrecondDiag(PreCond);
    jiba::rvec TomoCovar(TomoData.size());
    //we assume a general error of 5 ms for the seismic data
    std::fill(TomoCovar.begin(), TomoCovar.end(), 5.0);
    TomoObjective->SetDataCovar(TomoCovar);

    boost::shared_ptr<jiba::GravityObjective> ScalGravObjective(
        new jiba::GravityObjective(false, wantcuda));
    ScalGravObjective->SetObservedData(ScalGravData);
    ScalGravObjective->SetModelGeometry(GravModel);
    ScalGravObjective->SetPrecondDiag(PreCond);
    ScalGravObjective->SetDataCovar(jiba::ConstructError(ScalGravData, 0.02));

    boost::shared_ptr<jiba::GravityObjective> FTGObjective(
        new jiba::GravityObjective(true, wantcuda));
    FTGObjective->SetObservedData(FTGData);
    FTGObjective->SetModelGeometry(GravModel);
    FTGObjective->SetPrecondDiag(PreCond);
    FTGObjective->SetDataCovar(jiba::ConstructError(FTGData, 0.02));

    boost::shared_ptr<jiba::GradientRegularization> Regularization(
        new jiba::GradientRegularization(GravModel));

    Regularization->SetReferenceModel(InvModel);
    Regularization->SetDataCovar(InvModel);
    Regularization->SetPrecondDiag(PreCond);
    double tomolambda = 1.0;
    double scalgravlambda = 1.0;
    double ftglambda = 1.0;
    double reglambda = 1.0;
    std::cout << "Tomography Lambda: ";
    std::cin >> tomolambda;
    std::cout << "Scalar Gravimetry Lambda: ";
    std::cin >> scalgravlambda;
    std::cout << "FTG Lambda: ";
    std::cin >> ftglambda;
    std::cout << "Regularization Lambda: ";
    std::cin >> reglambda;
    size_t maxiter = 1;
    std::cout << "Maximum iterations: ";
    std::cin >> maxiter;

    boost::posix_time::ptime starttime =
        boost::posix_time::microsec_clock::local_time();
    if (tomolambda > 0.0)
      {
        Objective->AddObjective(TomoObjective, SlownessTransform, tomolambda);
      }
    if (scalgravlambda > 0.0)
      {
        Objective->AddObjective(ScalGravObjective, DensityTransform,
            scalgravlambda);
      }
    if (ftglambda > 0.0)
      {
        Objective->AddObjective(FTGObjective, DensityTransform, ftglambda);
      }
    Objective->AddObjective(Regularization, boost::shared_ptr<
        jiba::GeneralModelTransform>(new jiba::ModelCopyTransform()), reglambda);
    Objective->SetPrecondDiag(PreCond);
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

    jiba::rvec DensInvModel(DensityTransform->GeneralizedToPhysical(InvModel));
    jiba::rvec TomoInvModel(SlownessTransform->GeneralizedToPhysical(InvModel));

    std::copy(TomoInvModel.begin(), TomoInvModel.begin()
        + TomoModel.GetSlownesses().num_elements(),
        TomoModel.SetSlownesses().origin());
    std::copy(DensInvModel.begin(), DensInvModel.begin()
        + GravModel.SetDensities().num_elements(),
        GravModel.SetDensities().origin());

    //calculate the predicted refraction data
    std::cout << "Calculating response of inversion model." << std::endl;
    jiba::rvec TomoInvData(jiba::TomographyCalculator().Calculate(TomoModel));
    jiba::SaveTraveltimes(modelfilename + ".inv_tt.nc", TomoInvData, TomoModel);

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
    //and write out the data and model
    //here we have to distinguish again between scalar and ftg data
    std::cout << "Writing out inversion results." << std::endl;

    TomoModel.WriteNetCDF(modelfilename + ".tomo.inv.nc");
    TomoModel.WriteVTK(modelfilename + ".tomo.inv.vtk");
    GravModel.WriteVTK(modelfilename + ".grav.inv.vtk");
    GravModel.WriteNetCDF(modelfilename + ".grav.inv.nc");
    std::cout<<"Number of evaluations: " << Objective->GetNEval() << std::endl;
    boost::posix_time::ptime endtime =
        boost::posix_time::microsec_clock::local_time();
    double cachedruntime = (endtime - starttime).total_seconds();
    std::cout << "Runtime: " << cachedruntime << " s" << std::endl;
    std::cout << std::endl;

  }
