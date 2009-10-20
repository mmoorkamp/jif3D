//============================================================================
// Name        : tomogrid.cpp
// Author      : Max Moorkamp
// Version     :
// Copyright   : 2008, MM
//============================================================================


/*! \file tomogrid.cpp
 * A simple program to calculate evenly spaced refraction data from a model stored in a netcdf file.
 * The spacing and region for the measurements is specified interactively. The model is taken from
 * a netcdf file.
 */

#include "ThreeDSeismicModel.h"
#include "TomographyCalculator.h"
#include "ReadWriteTomographyData.h"
#include "../ModelBase/VTKTools.h"
#include "../Global/FileUtil.h"
#include "../Global/Noise.h"
#include <iostream>
#include <string>
#include <boost/cast.hpp>

using namespace std;

int main(int argc, char *argv[])
  {

    jiba::ThreeDSeismicModel SeisModel;

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
    //ask for the name of the netcdf file containing the model
    std::string ModelFilename = jiba::AskFilename("Model Filename: ");

    //determine the extension to find out the type
    std::string extension = jiba::GetFileExtension(ModelFilename);
    //read in the file
    if (extension == ".nc")
      {
        SeisModel.ReadNetCDF(ModelFilename);
      }
    else
      {
        std::cerr << "Wrong file format ! " << std::endl;
        return 100;
      }
    //setup the measurements in the forward modelling code
    const size_t nmeasx = boost::numeric_cast<size_t>((maxx - minx) / deltax);
    const size_t nmeasy = boost::numeric_cast<size_t>((maxy - miny) / deltay);
    for (size_t i = 0; i <= nmeasx; ++i)
      {
        for (size_t j = 0; j <= nmeasy; ++j)
          {
            SeisModel.AddMeasurementPoint(minx + i * deltax, miny + j * deltay,
                z);
            SeisModel.AddSource(minx + i * deltax, miny + j * deltay, z);
          }
      }
    const size_t nsource = SeisModel.GetSourcePosX().size();
    for (size_t i = 0; i < nsource; ++i)
      {
        for (size_t j = 0; j < nsource; ++j)
          {
            if (j != i)
              {
                SeisModel.AddMeasurementConfiguration(i, j);
              }
          }
      }

    jiba::TomographyCalculator Calculator;
    jiba::rvec TravelTimes(Calculator.Calculate(SeisModel));
    double error = 0.0;
    std::cout << "Traveltime error (ms): ";
    std::cin >> error;
    //if we want to add noise to the data
    if (error > 0.0)
      {
        jiba::AddNoise(TravelTimes, 0.0, error);
      }
    jiba::SaveTraveltimes(ModelFilename + ".tt.nc", TravelTimes, SeisModel);
    SeisModel.WriteVTK(ModelFilename + ".vtk");
  }
