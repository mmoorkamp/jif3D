//============================================================================
// Name        : wcompress.cpp
// Author      : Nov 7, 2008
// Version     :
// Copyright   : 2008, mmoorkamp
//============================================================================

#include <iostream>
#include <string>
#include <numeric>
#include "../Gravity/ScalarOMPGravityImp.h"
#include "../Gravity/ThreeDGravityModel.h"
#include "../Gravity/ReadWriteGravityData.h"
#include "../Global/Wavelet.h"
#include "../ModelBase/VTKTools.h"
#include "../Global/convert.h"
#include <boost/bind.hpp>
#include "../Gravity/WaveletCompressedGravityCalculator.h"

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
    //depending on the number of calling arguments
    switch (argc)
      {
    case 2:
      //1 argument, we assume measurement positions are store within the netcdf file for the model
      ModelFilename = argv[1];
      MeasPosFilename = argv[1];
      GravForward.ReadMeasPosNetCDF(MeasPosFilename);
      break;
    case 3:
      //2 arguments, we have a netcdf model and an ascii file with measurement positions
      ModelFilename = argv[1];
      MeasPosFilename = argv[2];
      GravForward.ReadMeasPosAscii(MeasPosFilename);
      break;
    default:
      //anything else, we ask for the filenames, measurement positions are ascii
      PromptForFiles(ModelFilename, MeasPosFilename);
      GravForward.ReadMeasPosAscii(MeasPosFilename);
      break;
      }
    //read the model from the netcdf file
    std::cout << "Reading model " << std::endl;
    GravForward.ReadNetCDF(ModelFilename);
    std::cout << "Calculating forward " << std::endl;

    boost::shared_ptr<jiba::WaveletCompressedGravityCalculator>
        WCalculator(jiba::CreateGravityCalculator<
            jiba::WaveletCompressedGravityCalculator>::MakeScalar());
    jiba::rvec Wavelet1(WCalculator->Calculate(GravForward));
    jiba::rvec Wavelet2(WCalculator->Calculate(GravForward));

    const size_t nmeas = GravForward.GetMeasPosX().size();
    for (size_t i = 0; i < nmeas; ++i)
      {

        std::cout << Wavelet1(i) << " " << Wavelet2(i) << " Rel. Error: "
            << (Wavelet1(i) - Wavelet2(i)) / Wavelet1(i) << std::endl;
      }
  }
