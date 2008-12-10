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
    boost::multi_array_types::size_type denssize[3];

    const size_t xsize = GravForward.GetDensities().shape()[0];
    const size_t ysize = GravForward.GetDensities().shape()[1];
    const size_t zsize = GravForward.GetDensities().shape()[2];
    const int nbglayers = GravForward.GetBackgroundDensities().size();

    denssize[0] = xsize;
    denssize[1] = ysize;
    denssize[2] = zsize; // + nbglayers;
    std::fill_n(transformsize, 3, 1.0);
    for (size_t i = 0; i < 3; ++i)
      {
        while (transformsize[i] < denssize[i])
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

    typedef boost::multi_array_types::index_range range;
    // OR typedef array_type::index_range range;

    size_t ncompressed = 0;

    for (size_t i = 0; i < nmeas; ++i)
      {
        boost::numeric::ublas::matrix_row<jiba::rmat> row(Sens, i);
        std::fill_n(CurrRow.origin(), CurrRow.num_elements(), 0.0);
        for (size_t j = 0; j < xsize; ++j)
          for (size_t k = 0; k < ysize; ++k)
            for (size_t l = 0; l < zsize; ++l)
              CurrRow[j][k][l] = row(j * (ysize * zsize) + k * zsize + l);

        //for (size_t j = 0; j < nbglayers; ++j)
        //  CurrRow[0][0][j] = row(xsize*ysize*zsize+j);

        jiba::WaveletTransform(CurrRow);
        size_t currn = 0;
        double relthresh = 1e-1;
        double normall = std::inner_product(CurrRow.origin(), CurrRow.origin()
            + CurrRow.num_elements(), CurrRow.origin(), 0.0);
        double normdiscarded = normall;
        const double maximum = *std::max_element(CurrRow.origin(),
            CurrRow.origin() + CurrRow.num_elements());
        while (normdiscarded / normall > 1e-2)
          {
            currn = 0;
            normdiscarded = 0.0;
            const double absthresh = maximum * relthresh;
            for (size_t j = 0; j < CurrRow.num_elements(); ++j)
              {
                if (fabs(*(CurrRow.origin() + j)) > absthresh)
                  {
                    SparseSens(i, j) = *(CurrRow.origin() + j);
                    currn++;
                  }
                else
                  {
                    normdiscarded += pow(*(CurrRow.origin() + j), 2);
                  }
              }
            relthresh /= 2.0;
            std::cout << "Measurement: " << i << " Threshold: " << relthresh
                << " Normratio: " << normdiscarded / normall << std::endl;
          }
        ncompressed += currn;
      }

    std::cout << "Compressing density vector " << std::endl;

    std::fill_n(CurrRow.origin(), CurrRow.num_elements(), 0.0);
    for (size_t j = 0; j < xsize; ++j)
      for (size_t k = 0; k < ysize; ++k)
        for (size_t l = 0; l < zsize; ++l)
          CurrRow[j][k][l] = GravForward.GetDensities()[j][k][l];

    //for (size_t j = 0; j < nbglayers; ++j)
    //  CurrRow[0][0][j] = GravForward.GetBackgroundDensities().at(j);
    jiba::WaveletTransform(CurrRow);

    jiba::rvec TransDens(CurrRow.num_elements());
    std::fill_n(TransDens.begin(), CurrRow.num_elements(), 0.0);
    std::copy(CurrRow.origin(), CurrRow.origin() + CurrRow.num_elements(),
        TransDens.begin());

    jiba::rvec SparseResult(boost::numeric::ublas::prec_prod(SparseSens,
        TransDens));

    std::cout << "Elements after compression: " << ncompressed << std::endl;
    //std::cout << SparseSens << std::endl;

    // calculate the size of the modelling domain for the background adjustment
    const double modelxwidth = std::accumulate(
        GravForward.GetXCellSizes().begin(), GravForward.GetXCellSizes().end(),
        0.0);
    const double modelywidth = std::accumulate(
        GravForward.GetYCellSizes().begin(), GravForward.GetYCellSizes().end(),
        0.0);
    const double modelzwidth = std::accumulate(
        GravForward.GetZCellSizes().begin(), GravForward.GetZCellSizes().end(),
        0.0);

    for (size_t i = 0; i < nmeas; ++i)
      {
        jiba::rmat Sens(xsize * ysize * zsize + nbglayers, 1);
        jiba::rvec BackResult(1);
        jiba::rmat dummy(0,0);

        BackResult = jiba::ScalarOMPGravityImp().CalcBackground(
            GravForward.GetMeasPosX()[i], GravForward.GetMeasPosY()[i],
            GravForward.GetMeasPosZ()[i], modelxwidth, modelywidth,
            modelzwidth, GravForward, dummy);
        SparseResult(i) += BackResult(0);
        std::cout << FullResult.at(i) << " " << SparseResult(i)
            << "Rel. Error: " << (FullResult.at(i) - SparseResult(i))
            / FullResult.at(i) << std::endl;
      }

    //GravForward.SaveScalarMeasurements(ModelFilename + ".out.nc");
    //GravForward.PlotScalarMeasurements(ModelFilename + ".plot");

    //GravForward.WriteVTK(ModelFilename + ".vtk");
  }
