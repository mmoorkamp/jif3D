//============================================================================
// Name        : wcompress.cpp
// Author      : Nov 7, 2008
// Version     :
// Copyright   : 2008, mmoorkamp
//============================================================================

#include <iostream>
#include <fstream>
#include <string>
#include <numeric>
#include "../Gravity/ScalarOMPGravityImp.h"
#include "../Gravity/ThreeDGravityModel.h"
#include "../Gravity/ReadWriteGravityData.h"
#include "../Gravity/DepthWeighting.h"
#include "../Gravity/ThreeDGravityFactory.h"
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

    jif3D::ThreeDGravityModel GravModel;
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
    std::ofstream outfile("accur.out");
    for (size_t i = 0; i < nruns; ++i)
      {
        boost::shared_ptr<jif3D::WaveletCompressedGravityCalculator>
            WCalculator(jif3D::CreateGravityCalculator<
                jif3D::WaveletCompressedGravityCalculator>::MakeTensor());

        const double z0 = 5.0;
        const double DepthExponent = -3.0;
        const size_t zsize = GravModel.GetDensities().shape()[2];
        jif3D::rvec WeightVector(zsize);

        jif3D::ConstructDepthWeighting(GravModel.GetZCellSizes(), z0,
            WeightVector, jif3D::WeightingTerm(DepthExponent));
        jif3D::rvec ModelWeight(GravModel.GetDensities().num_elements());
        ;
        std::fill_n(ModelWeight.begin(), ModelWeight.size(), 0);
        for (size_t i = 0; i < GravModel.GetDensities().num_elements(); ++i)
          {
            ModelWeight( i) = WeightVector(i % zsize);
          }
        WCalculator->SetWitheningVector().resize(ModelWeight.size());
        std::copy(ModelWeight.begin(),ModelWeight.end(),WCalculator->SetWitheningVector().begin());
        WCalculator->SetAccuracy(accuracy);
        jif3D::rvec Wavelet1(WCalculator->Calculate(GravModel));
        jif3D::rvec Wavelet2(WCalculator->Calculate(GravModel));

        const size_t nmeas = GravModel.GetMeasPosX().size();
        const size_t ndata = WCalculator->GetDataPerMeasurement();
        for (size_t i = 0; i < nmeas; ++i)
          {
            for (size_t j = 0; j < ndata; ++j)
              {
                const size_t index = i * ndata + j;
                std::cout << Wavelet1(index) << " " << Wavelet2(index)
                    << " Accuracy " << accuracy << " Rel. Error: "
                    << (Wavelet1(index) - Wavelet2(index)) / Wavelet1(index)
                    << std::endl;
              }
            std::cout << "Non zero: " << WCalculator->GetSensitivities().nnz()
                << std::endl;
            const double compression =
                double(WCalculator->GetSensitivities().nnz())
                    / double(GravModel.GetDensities().num_elements() * ndata
                        * nmeas);
            outfile << accuracy << " " << boost::numeric::ublas::norm_2(
                Wavelet1 - Wavelet2) / boost::numeric::ublas::norm_2(Wavelet1)
                << " " << compression << std::endl;
            std::cout << " Compression: " << compression << std::endl
                << std::endl;
          }
        accuracy /= 10.0;



      }

  }
