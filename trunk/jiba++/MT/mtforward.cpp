//============================================================================
// Name        : mtforward.cpp
// Author      : Jul 14, 2009
// Version     : 
// Copyright   : 2009, mmoorkamp
//============================================================================

#include <boost/numeric/conversion/cast.hpp>
#include "../Global/FileUtil.h"
#include "../Global/convert.h"
#include "../Global/Noise.h"
#include "X3DModel.h"
#include "X3DMTCalculator.h"
#include "ReadWriteImpedances.h"

int main()
  {
    std::string modelfilename = jiba::AskFilename("Filename: ", true);
    jiba::X3DModel MTModel;
    MTModel.ReadNetCDF(modelfilename);
    std::cout << "Model size: " << MTModel.GetXCellSizes().size() << " "
        << MTModel.GetYCellSizes().size() << " "
        << MTModel.GetZCellSizes().size() << std::endl;
    double minx, miny, maxx, maxy, deltax, deltay, z;
    //ask for the measurement grid specifications
    //first x-direction
    std::cout << "Minimum x-position: ";
    std::cin >> minx;
    std::cout << "Maximum x-position: ";
    std::cin >> maxx;
    std::cout << "Delta x: ";
    std::cin >> deltax;
    //then y-direction
    std::cout << "Minimum y-position: ";
    std::cin >> miny;
    std::cout << "Maximum y-position: ";
    std::cin >> maxy;
    std::cout << "Delta y: ";
    std::cin >> deltay;
    //all measurements have to be at the same height.
    std::cout << "Z-level: ";
    std::cin >> z;
    //setup the measurements in the forward modelling code
    const size_t nmeasx = boost::numeric_cast<size_t>((maxx - minx) / deltax);
    const size_t nmeasy = boost::numeric_cast<size_t>((maxy - miny) / deltay);
    for (size_t i = 0; i <= nmeasx; ++i)
      {
        for (size_t j = 0; j <= nmeasy; ++j)
          {
            MTModel.AddMeasurementPoint(minx + i * deltax, miny + j * deltay, z);
          }
      }

    std::string outfilename = jiba::AskFilename("Output filename: ", false);
    std::cout << "Frequencies: ";
    double currfreq = 1.0;
    std::vector<double> frequencies;
    try
      {
        while (true)
          {
            std::string input;
            std::cin >> input;
            jiba::convert(input, currfreq);
            frequencies.push_back(currfreq);
          }
      } catch (jiba::BadConversion &e)
      {

      }
    std::sort(frequencies.begin(), frequencies.end());
    std::copy(frequencies.begin(), frequencies.end(), std::back_inserter(
        MTModel.SetFrequencies()));
    std::cout << "Calculating forward response " << std::endl;
    jiba::X3DMTCalculator Calculator;
    jiba::rvec Impedances(Calculator.Calculate(MTModel));

    double relnoise = 0.0;
    std::cout << "Relative noise level: ";
    std::cin >> relnoise;
    jiba::AddNoise(Impedances, relnoise, 0.0);
    jiba::WriteImpedancesToNetCDF(outfilename, MTModel.GetFrequencies(),
        MTModel.GetMeasPosX(), MTModel.GetMeasPosY(), MTModel.GetMeasPosZ(),
        Impedances);
    MTModel.WriteVTK(modelfilename + ".vtk");
  }

