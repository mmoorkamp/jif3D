//============================================================================
// Name        : gravforward.cpp
// Author      : Max Moorkamp
// Version     :
// Copyright   : 2008, MM
//============================================================================


/*! \file gravforward.cpp
 * A simple program to calculate both scalar and tensor gravity data from a model stored in a netcdf file.
 * Measurement positions are also taken from a file. If one argument is given to the program,
 * the positions are taken from the netcdf file. If two arguments are given the first file is
 * a netcdf file containing the model, the second file an ascii file containing lines in the
 * form x,y,z that describe the position of the measurements in m. Otherwise the program will
 * ask for the name of a model file in netcdf format and the acii file described above.
 *
 * The program produces three output files with different endings, but with the model filename as a root. ".sdg.nc" is
 * a netcdf file with scalar gravity data. ".ftg.nc" is a netcdf file with tensor gravity data. ".vtk" is the model
 * in VTK format for plotting.
 */

#include <iostream>
#include <string>
#include "ThreeDGravityModel.h"
#include "ReadWriteGravityData.h"
#include "../GravMag/MinMemGravMagCalculator.h"
#include "ThreeDGravityFactory.h"
#include "../Global/FileUtil.h"



int main(int argc, char *argv[])
  {
    std::string ModelFilename, MeasPosFilename;

    jif3D::ThreeDGravityModel GravModel;
    //depending on the number of calling arguments
    switch (argc)
      {
    case 2:
      //1 argument, we assume measurement positions are stored within the netcdf file for the model
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
      ModelFilename = jif3D::AskFilename("Model Filename: ");
      MeasPosFilename = jif3D::AskFilename("Measurement Position Filename: ");
      GravModel.ReadMeasPosAscii(MeasPosFilename);
      break;
      }
    //read the model from the netcdf file
    GravModel.ReadNetCDF(ModelFilename);
    //create objects to calculate tensor and scalar data
    //we are not interested in sensitivities so we use a MinMem calculator to
    //allow for the calculation of large models.
    typedef typename jif3D::MinMemGravMagCalculator<jif3D::ThreeDGravityModel> CalculatorType;

    boost::shared_ptr<CalculatorType> TensorCalculator(jif3D::CreateGravityCalculator<CalculatorType>::MakeTensor());
    boost::shared_ptr<CalculatorType> ScalarCalculator(jif3D::CreateGravityCalculator<CalculatorType>::MakeScalar());

    //calculate both types of data
    jif3D::rvec ScalarResults(ScalarCalculator->Calculate(GravModel));
    jif3D::rvec TensorResults(TensorCalculator->Calculate(GravModel));
    jif3D::rvec ScalErr(ScalarResults.size(),0.0);
    jif3D::rvec TensErr(TensorResults.size(),0.0);
    //write the results to netcdf files
    jif3D::SaveScalarGravityMeasurements(ModelFilename + ".sgd.nc", ScalarResults, GravModel.GetMeasPosX(),
        GravModel.GetMeasPosY(), GravModel.GetMeasPosZ(), ScalErr);
    jif3D::SaveTensorGravityMeasurements(ModelFilename + ".ftg.nc", TensorResults, GravModel.GetMeasPosX(),
        GravModel.GetMeasPosY(), GravModel.GetMeasPosZ(),TensErr);
    //write the model in .vtk format, at the moment the best plotting option
    GravModel.WriteVTK(ModelFilename + ".vtk");
  }
