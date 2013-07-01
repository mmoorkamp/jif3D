//============================================================================
// Name        : gravinv.cpp
// Author      : Max Moorkamp
// Version     :
// Copyright   : 2008, MM
//============================================================================


/*! \file gravinv.cpp
 * Invert scalar or ftg gravity data. The program reads in a model file that specifies the starting model and
 * a file with the input data. It outputs the calculated data and the inversion model in netcdf and vtk file format
 * and the the raw and depth-weighted sensitivities in .vtk format for plotting.
 *
 * This program is a good example for the level of abstraction in jif3D++ and uses a lot of the current features. For the most
 * part it is completely transparent to the fact whether we invert scalar or FTG data. Only the reading and writing of data
 * requires knowledge of the data type.
 */

#include <iostream>
#include <fstream>
#include <string>
#include "../Global/convert.h"
#include "../Global/FatalException.h"
#include "../Global/FileUtil.h"
#include "../Inversion/LinearInversion.h"
#include "../ModelBase/VTKTools.h"
#include "../ModelBase/NetCDFModelTools.h"
#include "ThreeDGravityModel.h"
#include "ReadWriteGravityData.h"
#include "FullSensitivityGravityCalculator.h"
#include "DepthWeighting.h"
#include "ThreeDGravityFactory.h"

namespace ublas = boost::numeric::ublas;

//write the sensitivities for each measurement to a file
void WriteSensitivities(const std::string &nameroot,
    const std::string &sensname, const jif3D::rmat &Sens,
    const jif3D::ThreeDGravityModel &Model)
  {
    //create a data structure that mimics the geometry of the models
    jif3D::ThreeDModelBase::t3DModelData
        SensModel(
            boost::extents[Model.GetDensities().shape()[0]][Model.GetDensities().shape()[1]][Model.GetDensities().shape()[2]]);
    const size_t ngrid = Model.GetDensities().num_elements();
    const size_t ndata = Sens.size1();
    //for each measurement
    for (size_t i = 0; i < ndata; ++i)
      {
        //extract the corresponding row of the sensitivity matrix
        boost::numeric::ublas::matrix_row<const jif3D::rmat> sensrow(Sens, i);
        //copy to the Model structure to map to its geometric position
        std::copy(sensrow.begin(), sensrow.begin() + ngrid, SensModel.data());
        //write out a .vtk and a netcdf file
        jif3D::Write3DModelToVTK(nameroot + sensname + jif3D::stringify(i)
            + ".vtk", sensname, Model.GetXCellSizes(), Model.GetYCellSizes(),
            Model.GetZCellSizes(), SensModel);
        jif3D::Write3DModelToNetCDF(nameroot + sensname + jif3D::stringify(i)
            + ".nc", sensname, " ", Model.GetXCellSizes(),
            Model.GetYCellSizes(), Model.GetZCellSizes(), SensModel);
      }

  }

//calculate the normalized misfit between observed and synthetic data
double CalcMisfit(const jif3D::rvec &ObservedData,
    const jif3D::rvec &SyntheticData, jif3D::rvec &Misfit)
  {
    std::transform(ObservedData.begin(), ObservedData.end(),
        SyntheticData.begin(), Misfit.begin(), std::minus<double>());
    //we normalise the misfit by the observed data
    std::transform(Misfit.begin(), Misfit.end(), ObservedData.begin(),
        Misfit.begin(), std::divides<double>());
    return ublas::norm_2(Misfit);
  }

int main()
  {
    //these objects hold information about the measurements and their geometry
    jif3D::rvec Data;
    jif3D::ThreeDGravityModel::tMeasPosVec PosX, PosY, PosZ;

    //first we read in the starting model and the measured data
    std::string modelfilename = jif3D::AskFilename("Starting model Filename: ");
    std::string datafilename = jif3D::AskFilename( "Data Filename: ");

    //we read in the starting modelfile
    jif3D::ThreeDGravityModel Model;
    Model.ReadNetCDF(modelfilename);
    //get the name of the file containing the data and read it in

    //we figure out the type of data (scalar or ftg) from the variables
    //that are in the netcdf file
    jif3D::GravityDataType DataType = jif3D::IdentifyGravityDatafileType(
        datafilename);

    //Li and Oldenburg recommend a depth weighting exponent of -2.0
    //for scalar gravity data, we make this the default, it will be changed
    //for FTG
    double DepthExponent = -2.0;
    //create the pointer for the calculator object without assigning anything
    boost::shared_ptr<jif3D::FullSensitivityGravityCalculator> GravityCalculator;
    //now we have to do a few things differently depending on whether we deal
    //with scalar or FTG data
    //1. We have to read the data differently
    //2. We need a different forward modeling object
    //3. We need a different exponent for depth-weighting
    switch (DataType)
      {
    case jif3D::scalar:
      jif3D::ReadScalarGravityMeasurements(datafilename, Data, PosX, PosY, PosZ);
      //assign a scalar forward calculation object to the pointer
      GravityCalculator
          = boost::shared_ptr<jif3D::FullSensitivityGravityCalculator>(
              jif3D::CreateGravityCalculator<
                  jif3D::FullSensitivityGravityCalculator>::MakeScalar(true));
      break;
    case jif3D::ftg:
      jif3D::ReadTensorGravityMeasurements(datafilename, Data, PosX, PosY, PosZ);
      //assign a ftg forward calculation object to the pointer
      GravityCalculator
          = boost::shared_ptr<jif3D::FullSensitivityGravityCalculator>(
              jif3D::CreateGravityCalculator<
                  jif3D::FullSensitivityGravityCalculator>::MakeTensor(true));
      DepthExponent = -3.0;
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
    const size_t nmod = ngrid + Model.GetBackgroundDensities().size();

    //set the measurement points in the starting model to those of the data
    Model.ClearMeasurementPoints();
    for (size_t i = 0; i < nmeas; ++i)
      {
        Model.AddMeasurementPoint(PosX.at(i), PosY.at(i), PosZ.at(i));
      }
    //calculate the response of the starting model
    std::cout << "Calculating response of starting model." << std::endl;
    jif3D::rvec StartingData(GravityCalculator->Calculate(Model));

    std::cout << "Gridded model size: " << ngrid << " Complete size: " << nmod
        << std::endl;
    //we use this code to examine the behaviour of the sensitivities
    //so we write out the raw sensitivities for each measurement
    std::cout << "Writing out raw sensitivities." << std::endl;
    //WriteSensitivities(modelfilename, "raw_sens",
    //    GravityCalculator->GetSensitivities(), Model);

    //create objects for the depth weighting
    jif3D::rvec WeightVector(zsize), ModelWeight(nmod);
    //create objects for the misfit and a very basic error estimate
    jif3D::rvec DataDiffVec(ndata), DataError(ndata);

    //we create a simple error estimate by assuming 2% error
    //for each measurement
    std::cout << "Equalising sensitivity matrix." << std::endl;
    const double errorlevel = 0.02;
    double misfit = CalcMisfit(Data, StartingData, DataDiffVec);
    std::cout << "Initial misfit: " << misfit << std::endl;

    std::fill(DataError.begin(), DataError.end(), errorlevel);
    //and also equalise the sensitivity matrix
    for (size_t i = 0; i < ndata; ++i)
      {
        boost::numeric::ublas::matrix_row<jif3D::rmat> CurrentRow(
            GravityCalculator->SetSensitivities(), i);
        CurrentRow /= Data(i);
      }
    std::cout << "Calculating depth weighting." << std::endl;
    //now we perform the depth weighting for the sensitivities
    jif3D::rvec SensProfile;
    //we find a measurement site close to the centre of the model and extract the
    //sensitivity variation with depth
    jif3D::ExtractMiddleSens(Model, GravityCalculator->GetSensitivities(),
        GravityCalculator->GetDataPerMeasurement(), SensProfile);
    //we fit a curve of the form 1/(z+z0)^n to the extracted sensitivities
    double z0 = FitZ0(SensProfile, Model.GetZCellSizes(), jif3D::WeightingTerm(
        DepthExponent));
    std::cout << "Estimated z0: " << z0 << std::endl;

    //calculate the depth scaling
    jif3D::ConstructDepthWeighting(Model.GetZCellSizes(), z0, WeightVector,
        jif3D::WeightingTerm(DepthExponent));
    //and output the scaling weights
    //std::ofstream weightfile("weights.out");
    //std::copy(WeightVector.begin(), WeightVector.end(), std::ostream_iterator<
    //    double>(weightfile, "\n"));

    //we can choose the influence of the background on the inversion
    double backgroundweight = 1.0;
    std::cout << "Weight for background: ";
    std::cin >> backgroundweight;
    std::fill_n(ModelWeight.begin(), ModelWeight.size(), backgroundweight);

    //the WeightVector only has length zsize, one entry for each depth level
    //the inversion routine needs a vector with a weight for each model parameter
    // the weights only depend on the depth of the cell
    for (size_t i = 0; i < ngrid; ++i)
      {
        ModelWeight( i) = WeightVector(i % zsize);
      }

    //then we ask the user for the regularization parameter lambda
    double lambda = 1.0;
    std::cout << "Lambda: ";
    std::cin >> lambda;

    //here comes the core inversion
    //the problem is linear so we only perform a single inversion step
    std::cout << "Performing inversion." << std::endl;
    jif3D::rvec InvModel(nmod);
    std::fill(InvModel.begin(), InvModel.end(), 0.0);
    jif3D::DataSpaceInversion()(GravityCalculator->SetSensitivities(),
        DataDiffVec, ModelWeight, DataError, lambda, InvModel);

    //add the result of the inversion to the starting model
    //we only add the gridded part, the  background is always 0 due to the weighting
    std::transform(InvModel.begin(), InvModel.begin() + ngrid,
        Model.SetDensities().origin(), Model.SetDensities().origin(),
        std::plus<double>());
    jif3D::ThreeDGravityModel::tBackgroundVec Background(
        Model.GetBackgroundDensities());
    std::transform(InvModel.begin() + ngrid, InvModel.end(),
        Background.begin(), Background.begin(), std::plus<double>());
    Model.SetBackgroundDensities(Background);
    //calculate the predicted data
    std::cout << "Calculating response of inversion model." << std::endl;
    jif3D::rvec InvData(GravityCalculator->Calculate(Model));
    misfit = CalcMisfit(Data, InvData, DataDiffVec);
    std::cout << "Final misfit: " << misfit << std::endl;
    //and write out the data and model
    //here we have to distinguish again between scalar and ftg data
    std::cout << "Writing out inversion results." << std::endl;
    switch (DataType)
      {
    case jif3D::scalar:
      jif3D::Write3DDataToVTK(modelfilename + ".inv_sgd.vtk", "grav_accel",
          InvData, Model.GetMeasPosX(), Model.GetMeasPosY(),
          Model.GetMeasPosZ());
      jif3D::SaveScalarGravityMeasurements(modelfilename + ".inv_sgd.nc",
          InvData, Model.GetMeasPosX(), Model.GetMeasPosY(),
          Model.GetMeasPosZ());
      break;
    case jif3D::ftg:
      jif3D::Write3DTensorDataToVTK(modelfilename + ".inv_ftg.vtk",
          "grav_accel", InvData, Model.GetMeasPosX(), Model.GetMeasPosY(),
          Model.GetMeasPosZ());
      jif3D::SaveTensorGravityMeasurements(modelfilename + ".inv_ftg.nc",
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
    //WriteSensitivities(modelfilename, "fil_sens",
    //    GravityCalculator->GetSensitivities(), Model);
    std::cout << std::endl;
  }
