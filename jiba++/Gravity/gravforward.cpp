//============================================================================
// Name        : gravforward.cpp
// Author      : Max Moorkamp
// Version     :
// Copyright   : 2008, MM
//============================================================================


/*! \file A simple program to calculate gravity data from a model stored in a netcdf file.
 * Measurement positions are also taken from a file.
 * 
 */

#include <iostream>
#include <string>
#include "ThreeDGravityModel.h"

void PromptForFiles(std::string &ModelFilename, std::string &MeasPosFilename)
  {
    std::cout << "Model Filename: ";
    std::cin >> ModelFilename;
    std::cout << "Measurement Position Filename: ";
    std::cin >> MeasPosFilename;
  }

int main(int argc, char *argv[])
  {
    std::string ModelFilename, MeasPosFilename;
    
    jiba::ThreeDGravityModel GravForward;
    
    switch (argc)
      {
    case 2:
      ModelFilename = argv[1];
      MeasPosFilename = argv[1];
      GravForward.ReadMeasPosNetCDF(MeasPosFilename);
      break;
    case 3:
      ModelFilename = argv[1];
      MeasPosFilename = argv[2];
      GravForward.ReadMeasPosAscii(MeasPosFilename);
      break;
    default:
      PromptForFiles(ModelFilename, MeasPosFilename);
      GravForward.ReadMeasPosAscii(MeasPosFilename);
      break;
      }
    
    
    
    GravForward.ReadNetCDF(ModelFilename);
    GravForward.SaveScalarMeasurements(ModelFilename+".out.nc");
    GravForward.PlotScalarMeasurements(ModelFilename+".plot");
  }
