//============================================================================
// Name        : gravgrid.cpp
// Author      : Max Moorkamp
// Version     :
// Copyright   : 2008, MM
//============================================================================


/*! \file A simple program to calculate evenly spaced gravity data from a model stored in a netcdf file.
 * The spacing and region for the measurements is specified interactively. The model is taken from
 * a netcdf file.
 */

#include <iostream>
#include <string>
#include "ThreeDGravityModel.h"
#include "ReadWriteGravityData.h"
#include "../ModelBase/VTKTools.h"

using namespace std;

int main(int argc, char *argv[])
  {
    std::string ModelFilename;

    jiba::ThreeDGravityModel GravForward;

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
    const size_t nmeasx = (maxx - minx) / deltax;
    const size_t nmeasy = (maxy - miny) / deltay;
    for (size_t i = 0; i <= nmeasx; ++i)
      {
        for (size_t j = 0; j <= nmeasy; ++j)
          {
            GravForward.AddMeasurementPoint(minx + i * deltax, miny + j
                * deltay, z);
          }
      }
    //ask for the name of the netcdf file containing the model
    cout << "Model Filename: ";
    cin >> ModelFilename;
    //read in the file
    GravForward.ReadNetCDF(ModelFilename);
    //save the measurements and some plots
    jiba::ThreeDGravityModel::tScalarMeasVec Data(GravForward.CalcGravity());
    GravForward.CalcTensorGravity();
    GravForward.SaveScalarMeasurements(ModelFilename + ".out.nc");
    GravForward.SaveTensorMeasurements(ModelFilename + ".ftg.nc");
    GravForward.PlotScalarMeasurements(ModelFilename + ".plot");
    GravForward.PlotTensorMeasurements(ModelFilename);
    GravForward.WriteVTK(ModelFilename + ".vtk");
    jiba::Write3DDataToVTK(ModelFilename + ".data.vtk", "grav_accel", Data,
        GravForward.GetMeasPosX(), GravForward.GetMeasPosY(),
        GravForward.GetMeasPosZ());
  }
