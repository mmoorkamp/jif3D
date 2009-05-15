//============================================================================
// Name        : gravinv.cpp
// Author      : Max Moorkamp
// Version     :
// Copyright   : 2008, MM
//============================================================================


#include <iostream>
#include <fstream>
#include <string>
#include <boost/bind.hpp>
#include "../Global/convert.h"
#include "../Global/FatalException.h"
#include "../Global/NumUtil.h"
#include "../ModelBase/VTKTools.h"
#include "../ModelBase/NetCDFTools.h"
#include "../Gravity/ThreeDGravityModel.h"
#include "../Gravity/ReadWriteGravityData.h"
#include "../Gravity/FullSensitivityGravityCalculator.h"
#include "../Gravity/GravityObjective.h"
#include "../Gravity/DepthWeighting.h"
#include "../Gravity/GravityTransforms.h"
#include "../Inversion/LimitedMemoryQuasiNewton.h"
#include "../Inversion/NonLinearConjugateGradient.h"

namespace ublas = boost::numeric::ublas;

int main(int argc, char *argv[])
  {
    //these objects hold information about the measurements and their geometry
    jiba::rvec Data;
    jiba::ThreeDGravityModel::tMeasPosVec PosX, PosY, PosZ;

    //first we read in the starting model and the measured data
    std::string modelfilename, datafilename;
    std::cout << "Starting model Filename: ";
    std::cin >> modelfilename;
    //we read in the starting modelfile
    jiba::ThreeDGravityModel Model;
    Model.ReadNetCDF(modelfilename);
    //get the name of the file containing the data and read it in
    std::cout << "Data Filename: ";
    std::cin >> datafilename;

    //we figure out the type of data (scalar or ftg) from the variables
    //that are in the netcdf file
    jiba::GravityDataType DataType = jiba::IdentifyGravityDatafileType(
        datafilename);

    //Li and Oldenburg recommend a depth weighting exponent of -2.0
    //for scalar gravity data, we make this the default, it will be changed
    //for FTG
    double DepthExponent = -2.0;
    //create the pointer for the calculator object without assigning anything
    boost::shared_ptr<jiba::FullSensitivityGravityCalculator> GravityCalculator;
    boost::shared_ptr<jiba::VectorTransform> Transform;
    //now we have to do a few things differently depending on whether we deal
    //with scalar or FTG data
    //1. We have to read the data differently
    //2. We need a different forward modeling object
    //3. We need a different exponent for depth-weighting
    Transform = boost::shared_ptr<jiba::VectorTransform>(
        new jiba::CopyTransform());
    jiba::rvec InvarData;
    switch (DataType)
      {
    case jiba::scalar:
      jiba::ReadScalarGravityMeasurements(datafilename, Data, PosX, PosY, PosZ);
      //assign a scalar forward calculation object to the pointer
      GravityCalculator
          = boost::shared_ptr<jiba::FullSensitivityGravityCalculator>(
              jiba::CreateGravityCalculator<
                  jiba::FullSensitivityGravityCalculator>::MakeScalar());

      break;
    case jiba::ftg:
      jiba::ReadTensorGravityMeasurements(datafilename, Data, PosX, PosY, PosZ);
      //      //assign a ftg forward calculation object to the pointer
      GravityCalculator
          = boost::shared_ptr<jiba::FullSensitivityGravityCalculator>(
              jiba::CreateGravityCalculator<
                  jiba::FullSensitivityGravityCalculator>::MakeTensor());
      DepthExponent = -3.0;
      Transform = boost::shared_ptr<jiba::FTGInvariant>(
          new jiba::FTGInvariant());
      InvarData.resize((Data.size() / 9));
      for (size_t i = 0; i < Data.size(); i += 9)
        {
          jiba::rvec temp(Transform->Transform(ublas::vector_range<jiba::rvec>(
              Data, ublas::range(i, i + 9))));
          InvarData(i / 9) = temp(0);
        }
      Data.resize(InvarData.size());
      Data = InvarData;
      break;
    default:
      //in case we couldn't identify the data in the netcdf file
      //print an error message and exit with an error code
      std::cerr << "Cannot determine the type of data to invert. Aborting."
          << std::endl;
      exit(100);
      break;
      }
    //if we don't have data inversion doesn't make sense;
    if (Data.empty())
      {
        std::cerr << "No measurements defined" << std::endl;
        exit(100);
      }

    //we define a few constants that are used throughout the inversion
    const size_t nmeas = PosX.size();
    const size_t ndata = Data.size();
    const size_t xsize = Model.GetDensities().shape()[0];
    const size_t ysize = Model.GetDensities().shape()[1];
    const size_t zsize = Model.GetDensities().shape()[2];
    const size_t ngrid = xsize * ysize * zsize;
    const size_t nmod = ngrid;

    //set the measurement points in the starting model to those of the data
    Model.ClearMeasurementPoints();
    for (size_t i = 0; i < nmeas; ++i)
      {
        Model.AddMeasurementPoint(PosX.at(i), PosY.at(i), PosZ.at(i));
      }

    std::cout << "Gridded model size: " << ngrid << " Complete size: " << nmod
        << std::endl;

    //create objects for the misfit and a very basic error estimate
    jiba::rvec DataError(ndata);

    const double errorlevel = 0.02;
    const double maxdata = *std::max_element(Data.begin(), Data.end(),
        jiba::absLess<double, double>());
    for (size_t i = 0; i < ndata; ++i)
      {
        DataError( i) = std::max(std::abs(Data(i) * errorlevel), 1e-2 * maxdata
            * errorlevel);
      }
    std::cout << "DataError " << DataError << std::endl;

    boost::shared_ptr<jiba::GravityObjective> Objective(
        new jiba::GravityObjective(DataType == jiba::ftg));
    Objective->SetObservedData(Data);
    Objective->SetModelGeometry(Model);
    Objective->SetDataCovar(DataError);
    Objective->SetDataTransform(Transform);
    GravityCalculator->SetDataTransform(Transform);

    jiba::rvec InvModel(nmod);
    std::copy(Model.GetDensities().origin(), Model.GetDensities().origin()
        + Model.GetDensities().num_elements(), InvModel.begin());

    std::cout << "Calculating response of starting model." << std::endl;
    jiba::rvec StartingData(GravityCalculator->Calculate(Model));

    std::cout << "Calculating depth weighting." << std::endl;
    //now we perform the depth weighting for the sensitivities
    jiba::rvec SensProfile;
    jiba::rvec WeightVector(zsize), ModelWeight(nmod);
    //we find a measurement site close to the centre of the model and extract the
    //sensitivity variation with depth
    jiba::ExtractMiddleSens(Model, GravityCalculator->GetSensitivities(),
        GravityCalculator->GetDataPerMeasurement(), SensProfile);
    //we fit a curve of the form 1/(z+z0)^n to the extracted sensitivities
    double z0 = FitZ0(SensProfile, Model.GetZCellSizes(), jiba::WeightingTerm(
        DepthExponent));
    std::cout << "Estimated z0: " << z0 << std::endl;

    //calculate the depth scaling
    jiba::ConstructDepthWeighting(Model.GetZCellSizes(), z0, WeightVector,
        jiba::WeightingTerm(DepthExponent));
    for (size_t i = 0; i < ngrid; ++i)
      {
        ModelWeight( i) = WeightVector(i % zsize);
      }
    //here comes the core inversion

    std::cout << "Performing inversion." << std::endl;

    jiba::NonLinearConjugateGradient NLCG(Objective);
    jiba::LimitedMemoryQuasiNewton LBFGS(Objective, 5);

    LBFGS.SetModelCovDiag(ModelWeight);
    NLCG.SetModelCovDiag(ModelWeight);
    size_t iteration = 0;
    size_t maxiter = 30;
    do
      {
        LBFGS.MakeStep(InvModel);
        std::cout << std::endl;
        ++iteration;
      } while (iteration < maxiter && LBFGS.GetMisfit() > ndata);

    //add the result of the inversion to the starting model
    //we only add the gridded part, the  background is always 0 due to the weighting
    std::copy(InvModel.begin(), InvModel.begin() + ngrid,
        Model.SetDensities().origin());

    //calculate the predicted data
    std::cout << "Calculating response of inversion model." << std::endl;
    jiba::rvec InvData(GravityCalculator->Calculate(Model));
    std::cout << "Inversion Data: " << InvData << std::endl;
    //and write out the data and model
    //here we have to distinguish again between scalar and ftg data
    std::cout << "Writing out inversion results." << std::endl;
    switch (DataType)
      {
    case jiba::scalar:
      jiba::Write3DDataToVTK(modelfilename + ".inv_sgd.vtk", "grav_accel",
          InvData, Model.GetMeasPosX(), Model.GetMeasPosY(),
          Model.GetMeasPosZ());
      jiba::SaveScalarGravityMeasurements(modelfilename + ".inv_sgd.nc",
          InvData, Model.GetMeasPosX(), Model.GetMeasPosY(),
          Model.GetMeasPosZ());
      break;
    case jiba::ftg:
      //jiba::Write3DTensorDataToVTK(modelfilename + ".inv_ftg.vtk",
      //    "grav_accel", InvData, Model.GetMeasPosX(), Model.GetMeasPosY(),
      //    Model.GetMeasPosZ());
      //jiba::SaveTensorGravityMeasurements(modelfilename + ".inv_ftg.nc",
      //    InvData, Model.GetMeasPosX(), Model.GetMeasPosY(),
      //    Model.GetMeasPosZ());

      jiba::SaveScalarGravityMeasurements(modelfilename + ".meas_invariant.nc",
          Data, Model.GetMeasPosX(), Model.GetMeasPosY(), Model.GetMeasPosZ());
      jiba::SaveScalarGravityMeasurements(modelfilename + ".inv_invariant.nc",
          InvData, Model.GetMeasPosX(), Model.GetMeasPosY(),
          Model.GetMeasPosZ());
      break;
    default:
      std::cerr << " We should never reach this part. Fatal Error !"
          << std::endl;
      exit(100);
      break;
      }

    Model.WriteVTK(modelfilename + ".inv.vtk");
    Model.WriteNetCDF(modelfilename + ".inv.nc");
    std::cout << std::endl;
  }
