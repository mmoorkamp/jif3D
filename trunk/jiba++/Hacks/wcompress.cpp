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
    std::cout << "Reading model " << std::endl;
    GravModel.ReadNetCDF(ModelFilename);
    double accuracy = 0.1;
    const size_t nruns = 10;
    for (size_t i = 0; i < nruns; ++i)
      {
        boost::shared_ptr<jiba::WaveletCompressedGravityCalculator>
            WCalculator(jiba::CreateGravityCalculator<
                jiba::WaveletCompressedGravityCalculator>::MakeTensor());


        WCalculator->SetAccuracy(accuracy);
        jiba::rvec Wavelet1(WCalculator->Calculate(GravModel));
        jiba::rvec Wavelet2(WCalculator->Calculate(GravModel));

        const size_t nmeas = GravModel.GetMeasPosX().size();
        const size_t ndata = WCalculator->GetDataPerMeasurement();
        for (size_t i = 0; i < nmeas; ++i)
          {
            for (size_t j = 0; j < ndata; ++j)
              {
                const size_t index = i*ndata+j;
                std::cout << Wavelet1(index) << " " << Wavelet2(index) << " Accuracy "
                    << accuracy << " Rel. Error: " << (Wavelet1(index)
                    - Wavelet2(index)) / Wavelet1(index) << std::endl;
              }
            std::cout << " Compression: "
                << double(WCalculator->GetSensitivities().nnz())
                    / double(GravModel.GetDensities().num_elements())
                << std::endl << std::endl;
          }
        accuracy /= 10.0;
      }
  }
