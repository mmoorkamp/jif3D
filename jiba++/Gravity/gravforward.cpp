//============================================================================
// Name        : gravforward.cpp
// Author      : Max Moorkamp
// Version     :
// Copyright   : 2008, MM
//============================================================================


/*! \file gravforward.cpp
 * A simple program to calculate gravity data from a model stored in a netcdf file.
 * Measurement positions are also taken from a file. If one argument is given to the program,
 * the positions are taken from the netcdf file. If two arguments are given the first file is
 * a netcdf file containing the model, the second file an ascii file containing lines in the
 * form x,y,z that describe the position of the measurements in m. Otherwise the program will
 * ask for the name of a model file in netcdf format and the acii file described above.
 */

#include <iostream>
#include <string>
#include "ThreeDGravityModel.h"
#include "ReadWriteGravityData.h"
#include "MinMemGravityCalculator.h"

//get the name of the input files
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

    jiba::ThreeDGravityModel GravModel;
    //depending on the number of calling arguments
    switch (argc)
      {
    case 2:
      //1 argument, we assume measurement positions are store within the netcdf file for the model
      ModelFilename = argv[1];
      MeasPosFilename = argv[1];
      GravModel.ReadMeasPosNetCDF(MeasPosFilename);
      break;
    case 3:
      //2 arguments, we have a netcdf model and an ascii file with measurement positions
      ModelFilename = argv[1];
      MeasPosFilename = argv[2];
      GravModel.ReadMeasPosAscii(MeasPosFilename);
      break;
    default:
      //anything else, we ask for the filenames, measurement positions are ascii
      PromptForFiles(ModelFilename, MeasPosFilename);
      GravModel.ReadMeasPosAscii(MeasPosFilename);
      break;
      }
    //read the model from the netcdf file
    GravModel.ReadNetCDF(ModelFilename);
    boost::shared_ptr<jiba::MinMemGravityCalculator> TensorCalculator(jiba::CreateGravityCalculator<jiba::MinMemGravityCalculator>::MakeTensor());
    boost::shared_ptr<jiba::MinMemGravityCalculator> ScalarCalculator(jiba::CreateGravityCalculator<jiba::MinMemGravityCalculator>::MakeScalar());
    jiba::rvec ScalarResults(ScalarCalculator->Calculate(GravModel));
    jiba::rvec TensorResults(TensorCalculator->Calculate(GravModel));


    jiba::SaveScalarGravityMeasurements(ModelFilename + ".out.nc", ScalarResults, GravModel.GetMeasPosX(),
        GravModel.GetMeasPosY(), GravModel.GetMeasPosZ());
    jiba::SaveTensorGravityMeasurements(ModelFilename + ".ftg.nc", TensorResults, GravModel.GetMeasPosX(),
        GravModel.GetMeasPosY(), GravModel.GetMeasPosZ());
    //write the model in .vtk format, at the moment the best plotting option
    GravModel.WriteVTK(ModelFilename + ".vtk");
  }
