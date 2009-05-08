//============================================================================
// Name        : gravgrid.cpp
// Author      : Max Moorkamp
// Version     :
// Copyright   : 2008, MM
//============================================================================


/*! \file gravgrid.cpp
 * A simple program to calculate evenly spaced gravity data from a model stored in a netcdf file.
 * The spacing and region for the measurements is specified interactively. The model is taken from
 * a netcdf file.
 */

#include <iostream>
#include <string>
#include "ThreeDGravityModel.h"
#include "ReadWriteGravityData.h"
#include "MinMemGravityCalculator.h"
#include "../ModelBase/VTKTools.h"
#include "../Global/FileUtil.h"
#include <boost/cast.hpp>

using namespace std;

int main(int argc, char *argv[])
  {


    jiba::ThreeDGravityModel GravModel;

    double minx, miny, maxx, maxy, deltax, deltay, z;
    //ask for the measurement grid specifications
    //first x-direction
    cout << "Minimum x-position: ";
    cin >> minx;
    cout << "Maximum x-position: ";
    cin >> maxx;
    cout << "Delta x: ";
    cin >> deltax;
    //then y-direction
    cout << "Minimum y-position: ";
    cin >> miny;
    cout << "Maximum y-position: ";
    cin >> maxy;
    cout << "Delta y: ";
    cin >> deltay;
    //all measurements have to be at the same height.
    cout << "Z-level: ";
    cin >> z;
    //setup the measurements in the forward modelling code
    const size_t nmeasx = boost::numeric_cast<size_t>((maxx - minx) / deltax);
    const size_t nmeasy = boost::numeric_cast<size_t>((maxy - miny) / deltay);
    for (size_t i = 0; i <= nmeasx; ++i)
      {
        for (size_t j = 0; j <= nmeasy; ++j)
          {
            GravModel.AddMeasurementPoint(minx + i * deltax, miny + j
                * deltay, z);
          }
      }
    //ask for the name of the netcdf file containing the model
    std::string ModelFilename = jiba::AskFilename("Model Filename: ");

    //determine the extension to find out the type
    std::string extension = jiba::GetFileExtension(ModelFilename);
    //read in the file
    if (extension == ".nc")
      {
        GravModel.ReadNetCDF(ModelFilename);
      }
    else
      {
        GravModel.ReadIgmas(ModelFilename);
        GravModel.WriteNetCDF(ModelFilename + ".nc");
      }
    //save the measurements and some plots
    boost::shared_ptr<jiba::MinMemGravityCalculator>
        TensorCalculator(jiba::CreateGravityCalculator<
            jiba::MinMemGravityCalculator>::MakeTensor());
    boost::shared_ptr<jiba::MinMemGravityCalculator>
        ScalarCalculator(jiba::CreateGravityCalculator<
            jiba::MinMemGravityCalculator>::MakeScalar());
    jiba::rvec ScalarResults(ScalarCalculator->Calculate(GravModel));
    jiba::rvec TensorResults(TensorCalculator->Calculate(GravModel));

    jiba::SaveScalarGravityMeasurements(ModelFilename + ".sgd.nc",
        ScalarResults, GravModel.GetMeasPosX(), GravModel.GetMeasPosY(),
        GravModel.GetMeasPosZ());
    jiba::SaveTensorGravityMeasurements(ModelFilename + ".ftg.nc",
        TensorResults, GravModel.GetMeasPosX(), GravModel.GetMeasPosY(),
        GravModel.GetMeasPosZ());
    //write the model in .vtk format, at the moment the best plotting option
    GravModel.WriteVTK(ModelFilename + ".vtk");

    jiba::Write3DDataToVTK(ModelFilename + ".data.vtk", "grav_accel", ScalarResults,
        GravModel.GetMeasPosX(), GravModel.GetMeasPosY(),
        GravModel.GetMeasPosZ());
    jiba::Write3DTensorDataToVTK(ModelFilename + ".ftgdata.vtk", "U", TensorResults,
        GravModel.GetMeasPosX(), GravModel.GetMeasPosY(),
        GravModel.GetMeasPosZ());

  }
