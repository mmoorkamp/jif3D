//============================================================================
// Name        : jointinv.cpp
// Author      : May 12, 2009
// Version     :
// Copyright   : 2009, mmoorkamp
//============================================================================


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
#include <boost/bind.hpp>
#include "../Global/convert.h"
#include "../Global/FatalException.h"
#include "../Global/NumUtil.h"
#include "../Global/VectorTransform.h"
#include "../Global/FileUtil.h"
#include "../ModelBase/VTKTools.h"
#include "../ModelBase/NetCDFTools.h"
#include "../Inversion/LimitedMemoryQuasiNewton.h"
#include "../Inversion/JointObjective.h"
#include "../Inversion/MinDiffRegularization.h"
#include "../Tomo/ThreeDSeismicModel.h"
#include "../Tomo/ReadWriteTomographyData.h"
#include "../Tomo/TomographyObjective.h"
#include "../Tomo/TomographyCalculator.h"
#include "../Gravity/GravityObjective.h"
#include "../Gravity/ReadWriteGravityData.h"
#include "../Gravity/ThreeDGravityCalculator.h"
#include "../Gravity/MinMemGravityCalculator.h"
#include "../Gravity/DepthWeighting.h"

namespace ublas = boost::numeric::ublas;

class LogTransform: public jiba::GeneralModelTransform
  {
public:
  virtual jiba::rvec Transform(const jiba::rvec &FullModel)
    {
      jiba::rvec Output(FullModel.size());
      for (size_t i = 0; i < FullModel.size(); ++i)
        Output( i) = exp(FullModel(i));
      return Output;
    }

  virtual jiba::rvec Derivative(const jiba::rvec &FullModel,
      const jiba::rvec &Derivative)
    {

      jiba::rvec Output(FullModel.size());
      for (size_t i = 0; i < FullModel.size(); ++i)
        Output( i) = exp(FullModel(i)) * Derivative(i);

      return Output;
    }
  LogTransform()
    {
    }
  virtual ~LogTransform()
    {
    }
  };

class LogDensityTransform: public jiba::GeneralModelTransform
  {
public:
  virtual jiba::rvec Transform(const jiba::rvec &FullModel)
    {
      jiba::rvec Output(FullModel.size());
      for (size_t i = 0; i < FullModel.size(); ++i)
        Output( i) = (exp(-FullModel(i)) + 8500.0) / 5000.0;
      return Output;
    }

  virtual jiba::rvec Derivative(const jiba::rvec &FullModel,
      const jiba::rvec &Derivative)
    {
      const double factor = -1.0 / 5000.0;
      jiba::rvec Output(FullModel.size());
      for (size_t i = 0; i < FullModel.size(); ++i)
        {
          Output( i) = factor * exp(-FullModel(i)) * Derivative(i);
          // std::cout <<FullModel(i) << " " << Output(i) << " " << Derivative(i) << std::endl;
        }
      return Output;
    }
  LogDensityTransform()
    {
    }
  virtual ~LogDensityTransform()
    {
    }
  };

class VelTransform: public jiba::GeneralModelTransform
  {
  public:
  virtual jiba::rvec Transform(const jiba::rvec &FullModel)
    {
      jiba::rvec Output(FullModel.size());
      for (size_t i = 0; i < FullModel.size(); ++i)
        Output( i) = 1. / FullModel(i);
      return Output;
    }
  virtual jiba::rvec Derivative(const jiba::rvec &FullModel,
      const jiba::rvec &Derivative)
    {
      jiba::rvec Output(FullModel.size());
      for (size_t i = 0; i < FullModel.size(); ++i)
        Output( i) = -1. / (FullModel(i) * FullModel(i)) * Derivative(i);
      return Output;
    }
  VelTransform()
    {
    }
  virtual ~VelTransform()
    {
    }
  };

class VelDensTransform: public jiba::GeneralModelTransform
  {
  public:
  virtual jiba::rvec Transform(const jiba::rvec &FullModel)
    {
      jiba::rvec Output(FullModel.size());
      for (size_t i = 0; i < FullModel.size(); ++i)
        Output( i) = (FullModel(i) + 8500.0) / 5000.0;
      return Output;
    }
  virtual jiba::rvec Derivative(const jiba::rvec &FullModel,
      const jiba::rvec &Derivative)
    {
      jiba::rvec Output(FullModel.size());
      for (size_t i = 0; i < FullModel.size(); ++i)
        Output( i) = 5000.0 * Derivative(i);
      return Output;
    }
  VelDensTransform()
    {
    }
  virtual ~VelDensTransform()
    {
    }
  };

jiba::rvec ConstructError(const jiba::rvec &Data)
  {
    const size_t ndata = Data.size();
    const double errorlevel = 0.02;
    const double maxdata = *std::max_element(Data.begin(), Data.end(),
        jiba::absLess<double, double>());
    //create objects for the misfit and a very basic error estimate
    jiba::rvec DataError(ndata);
    for (size_t i = 0; i < ndata; ++i)
      {
        DataError( i) = std::max(std::abs(Data(i) * errorlevel), 1e-2 * maxdata
            * errorlevel);
      }
    return DataError;
  }

int main(int argc, char *argv[])
  {
    //these objects hold information about the measurements and their geometry
    jiba::rvec TomoData, GravData;

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

    std::string gravdatafilename = jiba::AskFilename("Gravity Data Filename: ");
    std::string gravmodelfilename = jiba::AskFilename(
        "Gravity Model Filename: ");
    jiba::ThreeDGravityModel GravModel;
    GravModel.ReadNetCDF(gravmodelfilename);
    GravModel = TomoModel;
    jiba::ThreeDGravityModel::tMeasPosVec PosX, PosY, PosZ;
    jiba::ReadScalarGravityMeasurements(gravdatafilename, GravData, PosX, PosY,
        PosZ);
    GravModel.ClearMeasurementPoints();
    for (size_t i = 0; i < PosX.size(); ++i)
      {
        GravModel.AddMeasurementPoint(PosX.at(i), PosY.at(i), PosZ.at(i));
      }

    //if we don't have data inversion doesn't make sense;
    if (TomoData.empty() || GravData.empty())
      {
        std::cerr << "No measurements defined" << std::endl;
        exit(100);
      }

    jiba::rvec InvModel(TomoModel.GetSlownesses().num_elements());
    std::copy(TomoModel.GetSlownesses().origin(),
        TomoModel.GetSlownesses().origin()
            + TomoModel.GetSlownesses().num_elements(), InvModel.begin());

    jiba::rvec RefModel(InvModel);
    for (size_t i = 0; i < InvModel.size(); ++i)
      InvModel( i) = 1./InvModel(i);
    boost::shared_ptr<jiba::TomographyObjective> TomoObjective(
        new jiba::TomographyObjective());
    TomoObjective->SetObservedData(TomoData);
    TomoObjective->SetModelGeometry(TomoModel);
    TomoObjective->SetDataCovar(ConstructError(TomoData));

    boost::shared_ptr<jiba::GravityObjective> GravObjective(
        new jiba::GravityObjective());
    GravObjective->SetObservedData(GravData);
    GravObjective->SetModelGeometry(GravModel);
    GravObjective->SetDataCovar(ConstructError(GravData));

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

    boost::shared_ptr<jiba::JointObjective> Objective(
        new jiba::JointObjective());
    boost::shared_ptr<jiba::MinDiffRegularization> Regularization(
        new jiba::MinDiffRegularization());

    Regularization->SetReferenceModel(RefModel);
    Regularization->SetDataCovar(RefModel);
    double gravlambda = 1.0;
    double reglambda = 1.0;
    std::cout << "Gravimetry Lambda: ";
    std::cin >> gravlambda;
    std::cout << "Regularization Lambda: ";
    std::cin >> reglambda;
    Objective->AddObjective(TomoObjective, boost::shared_ptr<
        jiba::GeneralModelTransform>(new VelTransform()));
    Objective->AddObjective(GravObjective, boost::shared_ptr<
        jiba::GeneralModelTransform>(new VelDensTransform()), gravlambda);
    Objective->AddObjective(Regularization, boost::shared_ptr<
        jiba::GeneralModelTransform>(new VelTransform()), reglambda);

    std::cout << "Performing inversion." << std::endl;

    jiba::LimitedMemoryQuasiNewton LBFGS(Objective, 5);
    LBFGS.SetModelCovDiag(ModelWeight);
    for (size_t i = 0; i < 5; ++i)
      {
        std::cout << "Iteration: " << i << std::endl;
        LBFGS.MakeStep(InvModel);
        TomoModel.WriteVTK(modelfilename + jiba::stringify(i) + ".tomo.inv.vtk");
        std::cout << std::endl;

      }

    jiba::rvec TomoInvModel(VelTransform().Transform(InvModel));
    jiba::rvec DensInvModel(VelDensTransform().Transform(InvModel));

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

    boost::shared_ptr<jiba::MinMemGravityCalculator>
        GravityCalculator =
            boost::shared_ptr<jiba::MinMemGravityCalculator>(
                jiba::CreateGravityCalculator<jiba::MinMemGravityCalculator>::MakeScalar(
                    true));
    jiba::rvec GravInvData(GravityCalculator->Calculate(GravModel));
    jiba::SaveScalarGravityMeasurements(modelfilename + ".inv_sgd.nc",
        GravInvData, GravModel.GetMeasPosX(), GravModel.GetMeasPosY(),
        GravModel.GetMeasPosZ());
    //and write out the data and model
    //here we have to distinguish again between scalar and ftg data
    std::cout << "Writing out inversion results." << std::endl;

    TomoModel.WriteVTK(modelfilename + ".tomo.inv.vtk");
    GravModel.WriteVTK(modelfilename + ".grav.inv.vtk");
    GravModel.WriteNetCDF(modelfilename + ".grav.inv.nc");
    std::cout << std::endl;
  }
