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

#include "../ModelBase/VTKTools.h"
#include "../Global/FileUtil.h"
#include "../Global/Noise.h"
#include <iostream>
#include <string>
#include <boost/cast.hpp>
#include "../Tomo/ReadWriteTomographyData.h"
#include "../Tomo/ThreeDSeismicModel.h"
#include "../Tomo/TomographyCalculator.h"
#include "../Tomo/TomographyData.h"

using namespace std;

int main(int argc, char *argv[])
  {

    jif3D::ThreeDSeismicModel SeisModel;
    jif3D::TomographyData Data;
    double recminx, recminy, recmaxx, recmaxy, recdeltax, recdeltay, recz;
    //ask for the measurement grid specifications
    //first x-direction
    std::cout << "Receiver configuration\n";
    cout << "  Minimum x-position: ";
    cin >> recminx;
    cout << "  Maximum x-position: ";
    cin >> recmaxx;
    cout << "  Delta x: ";
    cin >> recdeltax;
    //then y-direction
    cout << "  Minimum y-position: ";
    cin >> recminy;
    cout << "  Maximum y-position: ";
    cin >> recmaxy;
    cout << "  Delta y: ";
    cin >> recdeltay;
    //all measurements have to be at the same height.
    cout << "  Z-level: ";
    cin >> recz;

    double sorminx, sorminy, sormaxx, sormaxy, sordeltax, sordeltay, sorz;
    //ask for the source grid specifications
    //first x-direction
    std::cout << "Source configuration\n";
    cout << "  Minimum x-position: ";
    cin >> sorminx;
    cout << "  Maximum x-position: ";
    cin >> sormaxx;
    cout << "  Delta x: ";
    cin >> sordeltax;
    //then y-direction
    cout << "  Minimum y-position: ";
    cin >> sorminy;
    cout << "  Maximum y-position: ";
    cin >> sormaxy;
    cout << "  Delta y: ";
    cin >> sordeltay;
    //all measurements have to be at the same height.
    cout << "  Z-level: ";
    cin >> sorz;

    //ask for the name of the netcdf file containing the model
    std::string ModelFilename = jif3D::AskFilename("Model Filename: ");
    SeisModel.ReadNetCDF(ModelFilename);

    //setup the measurements in the forward modelling code
    const size_t recnmeasx = boost::numeric_cast<size_t>((recmaxx - recminx) / recdeltax);
    const size_t recnmeasy = boost::numeric_cast<size_t>((recmaxy - recminy) / recdeltay);
    for (size_t i = 0; i <= recnmeasx; ++i)
      {
        for (size_t j = 0; j <= recnmeasy; ++j)
          {
            Data.AddMeasurementPoint(recminx + i * recdeltax, recminy + j * recdeltay,
                recz);

          }
      }

    const size_t sornmeasx = boost::numeric_cast<size_t>((sormaxx - sorminx) / sordeltax);
    const size_t sornmeasy = boost::numeric_cast<size_t>((sormaxy - sorminy) / sordeltay);
    for (size_t i = 0; i <= sornmeasx; ++i)
      {
        for (size_t j = 0; j <= sornmeasy; ++j)
          {
            Data.AddSource(sorminx + i * sordeltax, sorminy + j * sordeltay, sorz);

          }
      }

    const size_t nsource = Data.GetSourcePosX().size();
    const size_t nrec = Data.GetMeasPosX().size();
    for (size_t i = 0; i < nsource; ++i)
      {
        for (size_t j = 0; j < nrec; ++j)
          {
            Data.AddMeasurementConfiguration(i, j);
          }
      }

    jif3D::TomographyCalculator Calculator;
    jif3D::rvec TravelTimes(Calculator.Calculate(SeisModel, Data));
    double error = 0.0;
    std::cout << "Traveltime error (s): ";
    std::cin >> error;
    //if we want to add noise to the data
    jif3D::rvec Errors(TravelTimes.size(), 0.0);
    if (error > 0.0)
      {
        jif3D::AddNoise(TravelTimes, 0.0, error);
        std::fill(Errors.begin(), Errors.end(), error);
      }

    Data.SetDataAndErrors(std::vector<double>(TravelTimes.begin(), TravelTimes.end()),
        std::vector<double>(Errors.begin(), Errors.end()));
    Data.WriteNetCDF(ModelFilename + ".tt.nc");

    SeisModel.WriteVTK(ModelFilename + ".vtk");
    Data.WriteMeasurementPoints(ModelFilename + ".rec.vtk");
    Data.WriteSourcePoints(ModelFilename + ".sor.vtk");

  }
