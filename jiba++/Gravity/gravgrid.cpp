//============================================================================
// Name        : gravgrid.cpp
// Author      : Max Moorkamp
// Version     :
// Copyright   : 2008, MM
//============================================================================


/*! \file A simple program to calculate evenly spaced gravity data from a model stored in a netcdf file.
 * The spacing and region for the measurements is specified interactively.
 * 
 */

#include <iostream>
#include <string>
#include "ThreeDGravityModel.h"

using namespace std;

int main(int argc, char *argv[])
  {
    std::string ModelFilename;

    jiba::ThreeDGravityModel GravForward;

    double minx, miny, maxx, maxy, deltax, deltay, z;

    cout << "Minimum x-position: ";
    cin >> minx;
    cout << "Maximum x-position: ";
    cin >> maxx;
    cout << "Delta x: ";
    cin >> deltax;

    cout << "Minimum y-position: ";
    cin >> miny;
    cout << "Maximum y-position: ";
    cin >> maxy;
    cout << "Delta y: ";
    cin >> deltay;
    
    cout << "Z-level: ";
    cin >> z;
    
    cout << "Model Filename: ";
    cin >> ModelFilename;
    const size_t nmeasx = (maxx - minx)/deltax;
    const size_t nmeasy = (maxy - miny)/deltay;
    for (size_t i = 0; i < nmeasx; ++i)
      {
        for (size_t j = 0; j < nmeasy; ++j) {
                GravForward.AddMeasurementPoint(minx+i*deltax,miny+j*deltay,z);
        }
      }
    GravForward.ReadNetCDF(ModelFilename);
    GravForward.SaveScalarMeasurements(ModelFilename+".out.nc");
    GravForward.PlotScalarMeasurements(ModelFilename+".plot");
  }
