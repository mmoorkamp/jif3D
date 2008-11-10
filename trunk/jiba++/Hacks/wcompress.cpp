//============================================================================
// Name        : wcompress.cpp
// Author      : Nov 7, 2008
// Version     :
// Copyright   : 2008, mmoorkamp
//============================================================================

#include <iostream>
#include <string>
#include <numeric>
#include "../Gravity/ThreeDGravityModel.h"
#include "../Gravity/ReadWriteGravityData.h"
#include "../Global/Wavelet.h"
#include <boost/bind.hpp>

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

    jiba::ThreeDGravityModel GravForward(true);
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
    jiba::ThreeDGravityModel::tScalarMeasVec FullResult(
        GravForward.CalcGravity());
    std::cout << "Getting sensitivities " << std::endl;
    jiba::rmat Sens(GravForward.GetScalarSensitivities());

    std::cout << "Compressing matrix " << std::endl;
    const size_t nmeas = Sens.size1();
    const size_t nparam = Sens.size2();
    boost::multi_array_types::size_type transformsize[3];
    std::fill_n(transformsize, 3, 1.0);
    for (size_t i = 0; i < 3; ++i)
      {
        while (transformsize[i] < GravForward.GetDensities().shape()[i])
          {
            transformsize[i] *= 2;
          }
        std::cout << "Transformsize Dim: " << i << " Size: "
            << transformsize[i] << std::endl;
      }
    boost::multi_array<double, 3> CurrRow(
        boost::extents[transformsize[0]][transformsize[1]][transformsize[2]]);
    boost::numeric::ublas::mapped_matrix<double,
        boost::numeric::ublas::column_major> SparseSens(nmeas,
        CurrRow.num_elements());
    const double relthresh = 1e-3;
    size_t ncompressed = 0;
    for (size_t i = 0; i < nmeas; ++i)
      {
        boost::numeric::ublas::matrix_row<jiba::rmat> row(Sens, i);
        std::fill_n(CurrRow.origin(), CurrRow.num_elements(), 0.0);
        copy(row.begin(), row.end(), CurrRow.origin());
        jiba::WaveletTransform(CurrRow);
        const double absthresh = *std::max_element(CurrRow.origin(),
            CurrRow.origin() + CurrRow.num_elements()) * relthresh;

        for (size_t j = 0; j < CurrRow.num_elements(); ++j)
          {
            if (fabs(*(CurrRow.origin() + j)) > absthresh)
              {
                SparseSens(i, j) = *(CurrRow.origin() + j);
                ncompressed++;
              }
          }
      }

    const size_t xsize = GravForward.GetDensities().shape()[0];
    const size_t ysize = GravForward.GetDensities().shape()[1];
    const size_t zsize = GravForward.GetDensities().shape()[2];



    std::cout << "Compressing density vector " << std::endl;

    std::fill_n(CurrRow.origin(), CurrRow.num_elements(), 0.0);
    std::copy(GravForward.GetDensities().origin(),
        GravForward.GetDensities().origin()
            + GravForward.GetDensities().num_elements(), CurrRow.origin());
    copy(GravForward.GetBackgroundDensities().begin(),
        GravForward.GetBackgroundDensities().end(), CurrRow.origin() + xsize
            * ysize * zsize);
    jiba::WaveletTransform(CurrRow);

    jiba::rvec TransDens(CurrRow.num_elements());
    std::fill_n(TransDens.begin(), CurrRow.num_elements(), 0.0);
    std::copy(CurrRow.origin(), CurrRow.origin() + CurrRow.num_elements(),
        TransDens.begin());

    jiba::rvec SparseResult(boost::numeric::ublas::prec_prod(SparseSens,
        TransDens));

    std::cout << "Elements after compression: " << ncompressed << std::endl;
    //std::cout << SparseSens << std::endl;

    for (size_t i = 0; i < nmeas; ++i)
      std::cout << FullResult.at(i) << " " << SparseResult(i) << std::endl;

    //GravForward.SaveScalarMeasurements(ModelFilename + ".out.nc");
    //GravForward.PlotScalarMeasurements(ModelFilename + ".plot");

    //GravForward.WriteVTK(ModelFilename + ".vtk");
  }
