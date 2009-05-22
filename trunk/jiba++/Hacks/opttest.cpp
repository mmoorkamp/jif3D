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
#include "OptLBFGS.h"
#include "NLF.h"
#include "../Gravity/MinMemGravityCalculator.h"
#include "../Inversion/JointObjective.h"
#include <boost/bind.hpp>
#include "../Global/convert.h"
#include "../Global/FatalException.h"
#include "../Global/NumUtil.h"
#include "../Global/VectorTransform.h"
#include "../Global/FileUtil.h"
#include "../Inversion/ModelTransforms.h"
#include "../Inversion/JointObjective.h"
#include "../Inversion/MinDiffRegularization.h"
#include "../Inversion/ConstructError.h"
#include "../ModelBase/VTKTools.h"
#include "../ModelBase/NetCDFTools.h"
#include "../Gravity/GravityObjective.h"
#include "../Gravity/ReadWriteGravityData.h"
#include "../Gravity/ThreeDGravityCalculator.h"
#include "../Gravity/MinMemGravityCalculator.h"
#include "../Gravity/DepthWeighting.h"
using NEWMAT::ColumnVector;

boost::shared_ptr<jiba::JointObjective> Objective(new jiba::JointObjective());
jiba::rvec InvModel;

void update_model(int, int, ColumnVector x)
  {
    for (size_t i = 0; i < InvModel.size(); ++i)
      InvModel(i) = x(i+1);
  }

void init_model(int ndim, ColumnVector& x)
  {
    if (ndim != InvModel.size())
      {
        exit(1);
      }
    for (size_t i = 0; i < InvModel.size(); ++i)
      x(i + 1) = InvModel(i);
  }

void eval_objective(int mode, int n, const ColumnVector& x, double& fx,
    ColumnVector& g, int& result)
  {

    if (n != InvModel.size())
      return;

    for (size_t i = 0; i < InvModel.size(); ++i)
      InvModel( i) = x(i + 1);

    if (mode & OPTPP::NLPFunction)
      {
        fx = Objective->CalcMisfit(InvModel);
        result = OPTPP::NLPFunction;
      }
    if (mode & OPTPP::NLPGradient)
      {
        jiba::rvec Gradient(Objective->CalcGradient());
        for (size_t i = 0; i < Gradient.size(); ++i)
          g(i + 1) = Gradient(i);

        result = OPTPP::NLPGradient;
      }
  }



int main(int argc, char *argv[])
  {
    //these objects hold information about the measurements and their geometry
    jiba::rvec ScalGravData, FTGData;
    jiba::ThreeDGravityModel GravModel;
    //first we read in the starting model and the measured data
    std::string modelfilename = jiba::AskFilename("Starting model Filename: ");
    GravModel.ReadNetCDF(modelfilename);

    std::string scalgravdatafilename = jiba::AskFilename(
        "Scalar Gravity Data Filename: ");
    std::string ftgdatafilename = jiba::AskFilename("FTG Data Filename: ");

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
    if (ScalGravData.empty())
      {
        std::cerr << "No measurements defined" << std::endl;
        exit(100);
      }

    InvModel.resize(GravModel.GetDensities().num_elements());
    std::copy(GravModel.GetDensities().origin(),
        GravModel.GetDensities().origin()
            + GravModel.GetDensities().num_elements(), InvModel.begin());

    jiba::rvec RefModel(InvModel);
    std::fill(RefModel.begin(), RefModel.end(), 1.0);
    boost::shared_ptr<jiba::GeneralModelTransform> Transform(
        new jiba::ModelCopyTransform);
    InvModel = Transform->PhysicalToGeneralized(InvModel);

    boost::shared_ptr<jiba::GravityObjective> ScalGravObjective(
        new jiba::GravityObjective());
    ScalGravObjective->SetObservedData(ScalGravData);
    ScalGravObjective->SetModelGeometry(GravModel);
    ScalGravObjective->SetDataCovar(jiba::ConstructError(ScalGravData, 0.02));

    boost::shared_ptr<jiba::GravityObjective> FTGObjective(
        new jiba::GravityObjective(true));
    FTGObjective->SetObservedData(FTGData);
    FTGObjective->SetModelGeometry(GravModel);
    FTGObjective->SetDataCovar(jiba::ConstructError(FTGData, sqrt(0.02)));

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

    boost::shared_ptr<jiba::MinDiffRegularization> Regularization(
        new jiba::MinDiffRegularization());

    Regularization->SetReferenceModel(RefModel);
    Regularization->SetDataCovar(RefModel);

    double scalgravlambda = 1.0;
    double ftglambda = 1.0;
    double reglambda = 1.0;

    std::cout << "Scalar Gravimetry Lambda: ";
    std::cin >> scalgravlambda;
    std::cout << "FTG Lambda: ";
    std::cin >> ftglambda;
    std::cout << "Regularization Lambda: ";
    std::cin >> reglambda;

    if (scalgravlambda > 0.0)
      {
        Objective->AddObjective(ScalGravObjective, Transform, scalgravlambda);
      }
    if (ftglambda > 0.0)
      {
        Objective->AddObjective(FTGObjective, Transform, ftglambda);
      }
    Objective->AddObjective(Regularization, Transform, reglambda);

    int n = InvModel.size();

    static char *status_file =
      { "tstLBFGS.out" };

    //  Create a Nonlinear problem object

    OPTPP::NLF1 nlp(n, eval_objective, init_model);

    //  Build a LBFGS object and optimize

    OPTPP::OptLBFGS objfcn(&nlp);
    objfcn.setUpdateModel(update_model);
    if (!objfcn.setOutputFile(status_file, 0))
      cerr << "main: output file open failed" << endl;
    objfcn.setGradTol(1.e-6);
    objfcn.setMaxBacktrackIter(10);
    objfcn.setPrintFinalX(true);
    objfcn.optimize();

    objfcn.printStatus("Solution from LBFGS: More and Thuente's linesearch");

    objfcn.cleanup();

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
    //and write out the data and model
    //here we have to distinguish again between scalar and ftg data
    std::cout << "Writing out inversion results." << std::endl;
    GravModel.WriteVTK(modelfilename + ".grav.inv.vtk");
    GravModel.WriteNetCDF(modelfilename + ".grav.inv.nc");
    std::cout << std::endl;

  }
