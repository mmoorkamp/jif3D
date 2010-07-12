//============================================================================
// Name        : gravconst.cpp
// Author      : Max Moorkamp
// Version     :
// Copyright   : 2008, MM
//============================================================================


#include <iostream>
#include <fstream>
#include <string>
#include "../Global/convert.h"
#include "../Global/FatalException.h"
#include "../Global/FileUtil.h"
#include "../Global/NumUtil.h"
#include "../Inversion/LinearInversion.h"
#include "../Inversion/MatrixTools.h"
#include "../ModelBase/VTKTools.h"
#include "../ModelBase/NetCDFTools.h"
#include "ThreeDGravityModel.h"
#include "ReadWriteGravityData.h"
#include "FullSensitivityGravityCalculator.h"
#include "DepthWeighting.h"
#include "ThreeDGravityFactory.h"

namespace ublas = boost::numeric::ublas;

//calculate the normalized misfit between observed and synthetic data
double CalcMisfit(const jiba::rvec &ObservedData,
    const jiba::rvec &SyntheticData, jiba::rvec &Misfit)
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
    jiba::rvec Data;
    jiba::ThreeDGravityModel::tMeasPosVec PosX, PosY, PosZ;

    std::string datafilename = jiba::AskFilename("Data Filename: ");

    //get the name of the file containing the data and read it in

    //we figure out the type of data (scalar or ftg) from the variables
    //that are in the netcdf file
    jiba::GravityDataType DataType = jiba::IdentifyGravityDatafileType(
        datafilename);

    //create the pointer for the calculator object without assigning anything
    boost::shared_ptr<jiba::FullSensitivityGravityCalculator> GravityCalculator;
    //now we have to do a few things differently depending on whether we deal
    //with scalar or FTG data
    //1. We have to read the data differently
    //2. We need a different forward modeling object
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
      //assign a ftg forward calculation object to the pointer
      GravityCalculator
          = boost::shared_ptr<jiba::FullSensitivityGravityCalculator>(
              jiba::CreateGravityCalculator<
                  jiba::FullSensitivityGravityCalculator>::MakeTensor());
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

    const double maxxmeas = *std::max_element(PosX.begin(), PosX.end());
    const double maxymeas = *std::max_element(PosY.begin(), PosY.end());
    const double minxmeas = *std::min_element(PosX.begin(), PosX.end());
    const double minymeas = *std::min_element(PosY.begin(), PosY.end());

    jiba::rvec XPosDiff(nmeas), YPosDiff(nmeas);
    std::adjacent_difference(PosX.begin(), PosX.end(), XPosDiff.begin());
    std::adjacent_difference(PosY.begin(), PosY.end(), YPosDiff.begin());

    const int nrefine = 3;
    const double deltax = (maxxmeas - minxmeas)/(nrefine * sqrt(double(ndata)));
    const double deltay = (maxymeas - minymeas)/(nrefine * sqrt(double(ndata)));

    std::cout << "Delta X: " << deltax << " Delta Y: " << deltay << std::endl;
    const int nx = std::ceil((maxxmeas - minxmeas) / deltax);
    const int ny = std::ceil((maxymeas  - minymeas)/ deltay);
    const int nz = 1;
    const double deltaz = 1.0;
    std::cout << "Nx: " << nx << " Ny: " << ny << std::endl;
    jiba::ThreeDGravityModel Model;

    Model.SetXCellSizes().resize(boost::extents[nx]);
    Model.SetYCellSizes().resize(boost::extents[ny]);
    Model.SetZCellSizes().resize(boost::extents[nz]);
    Model.SetDensities().resize(boost::extents[nx][ny][nz]);
    double defaultdensity = 0.0;
    std::fill_n(Model.SetDensities().origin(),
        Model.GetDensities().num_elements(), defaultdensity);
    std::fill_n(Model.SetXCellSizes().begin(), nx, deltax);
    std::fill_n(Model.SetYCellSizes().begin(), ny, deltay);
    std::fill_n(Model.SetZCellSizes().begin(), nz, deltaz);
    Model.SetOrigin(minxmeas, minymeas, 100.0);

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
    //we need this to get the sensitivities
    std::cout << "Calculating response of starting model." << std::endl;
    jiba::rvec StartingData(GravityCalculator->Calculate(Model));
    std::cout << "Gridded model size: " << ngrid << " Complete size: " << nmod
        << std::endl;

    //create objects for the depth weighting
    jiba::rvec ModelWeight(nmod);
    std::fill(ModelWeight.begin(), ModelWeight.end(), 1.0);
    //create objects for the misfit and a very basic error estimate
    jiba::rvec DataDiffVec(ndata), DataError(ndata);

    //we create a simple error estimate by assuming 2% error
    //for each measurement
    std::cout << "Equalizing sensitivity matrix." << std::endl;
    const double errorlevel = 0.02;
    double misfit = CalcMisfit(Data, StartingData, DataDiffVec);
    std::cout << "Initial misfit: " << misfit << std::endl;

    std::fill(DataError.begin(), DataError.end(), errorlevel);
    //and also equalise the sensitivity matrix
    for (size_t i = 0; i < ndata; ++i)
      {
        boost::numeric::ublas::matrix_row<jiba::rmat> CurrentRow(
            GravityCalculator->SetSensitivities(), i);
        CurrentRow /= Data(i);
      }

    //here comes the core inversion
    //the problem is linear so we only perform a single inversion step
    std::cout << "Performing inversion." << std::endl;
    jiba::rvec InvModel(nmod);
    std::fill(InvModel.begin(), InvModel.end(), 0.0);
    jiba::rmat InvMat;
    jiba::GeneralizedInverse()(GravityCalculator->SetSensitivities(), InvMat,
        1e-14);
    jiba::rvec Ones(ndata);
    std::fill_n(Ones.begin(), ndata, 1.0);
    InvModel = boost::numeric::ublas::prod(InvMat, Ones);
    std::copy(InvModel.begin(), InvModel.end(), Model.SetDensities().origin());

    //calculate the predicted data
    double newgrid = 0.0, newz = 0.0;
    std::cout << "New grid spacing: ";
    std::cin >> newgrid;
    std::cout << "New z-level: ";
    std::cin >> newz;
    Model.ClearMeasurementPoints();
    for (double x = minxmeas; x < maxxmeas; x += newgrid)
      {
        for (double y = minymeas; y < maxymeas; y += newgrid)
          Model.AddMeasurementPoint(x, y, newz);
      }
    std::cout << "Calculating response of inversion model." << std::endl;
    jiba::rvec InvData(GravityCalculator->Calculate(Model));
    //and write out the data and model
    //here we have to distinguish again between scalar and ftg data
    std::string modelfilename = "pseudo";
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
      jiba::Write3DTensorDataToVTK(modelfilename + ".inv_ftg.vtk",
          "grav_accel", InvData, Model.GetMeasPosX(), Model.GetMeasPosY(),
          Model.GetMeasPosZ());
      jiba::SaveTensorGravityMeasurements(modelfilename + ".inv_ftg.nc",
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
