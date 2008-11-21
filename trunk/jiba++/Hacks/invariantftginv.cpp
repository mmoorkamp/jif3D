//============================================================================
// Name        : ftginv.cpp
// Author      : Max Moorkamp
// Version     :
// Copyright   : 2008, MM
//============================================================================


#include <iostream>
#include <fstream>
#include <string>
#include "../Global/convert.h"
#include "ThreeDGravityModel.h"
#include "../Inversion/MatrixTools.h"
#include "ReadWriteGravityData.h"
#include <boost/numeric/bindings/atlas/cblas2.hpp>
#include "../Inversion/LinearInversion.h"
#include "../ModelBase/VTKTools.h"
#include "DepthWeighting.h"

namespace atlas = boost::numeric::bindings::atlas;
namespace ublas = boost::numeric::ublas;

double CalcInvariant(const jiba::GravimetryMatrix &Matrix)
  {
    return Matrix(0, 0) * Matrix(1, 1) + Matrix(1, 1) * Matrix(2, 2) + Matrix(
        0, 0) * Matrix(2, 2) - Matrix(0, 1) * Matrix(0, 1) - Matrix(1, 2)
        * Matrix(1, 2) - Matrix(0, 2) * Matrix(0, 2);
  }

int main(int argc, char *argv[])
  {
    jiba::ThreeDGravityModel Model(false, true);

    jiba::ThreeDGravityModel::tTensorMeasVec Data;
    jiba::ThreeDGravityModel::tMeasPosVec PosX, PosY, PosZ;

    std::string modelfilename, datafilename;
    std::cout << "Mesh Filename: ";
    std::cin >> modelfilename;
    //we read in a complete modelfile, but we only use the mesh information
    Model.ReadNetCDF(modelfilename);
    //get the name of the file containing the data and read it in
    std::cout << "Data Filename: ";
    std::cin >> datafilename;
    jiba::ReadTensorGravityMeasurements(datafilename, Data, PosX, PosY, PosZ);
    double lambda = 1.0;
    std::cout << "Lambda: ";
    std::cin >> lambda;
    const size_t ndata = Data.size();
    jiba::rvec DataVector(ndata), DataError(ndata);
    for (size_t i = 0; i < ndata; ++i)
      {
        DataVector( i) = CalcInvariant(Data.at(i));
      }

    const size_t nmeas = PosX.size();
    const size_t xsize = Model.GetDensities().shape()[0];
    const size_t ysize = Model.GetDensities().shape()[1];
    const size_t zsize = Model.GetDensities().shape()[2];
    const size_t ngrid = xsize * ysize * zsize;
    const size_t nmod = ngrid + Model.GetBackgroundDensities().size();
    jiba::rvec InvModel(nmod), StartingModel(nmod);
    std::fill_n(InvModel.begin(), nmod, 0.0);
    std::copy(Model.GetDensities().origin(), Model.GetDensities().origin()
        + ngrid, StartingModel.begin());
    std::copy(Model.GetBackgroundDensities().begin(),
        Model.GetBackgroundDensities().end(), StartingModel.begin() + ngrid);
    //set the measurement points in the model to those of the data
    Model.ClearMeasurementPoints();
    for (size_t i = 0; i < nmeas; ++i)
      {
        Model.AddMeasurementPoint(PosX.at(i), PosY.at(i), PosZ.at(i));
      }

    const double errorlevel = 0.02;
    std::transform(DataVector.begin(), DataVector.end(), DataError.begin(),
        boost::bind(std::abs<double>, boost::bind(std::multiplies<double>(),
            _1, errorlevel)));
    for (size_t i = 0; i < DataError.size(); ++i)
      {
        DataError( i) = std::max(DataError(i), 1e-25);
      }

    const size_t maxiter = 30;
    size_t iter = 0;
    double gradnorm = 1e10;
    const double mingradnorm = 0.1;
    double misfit = 1e10;
    const double minmisfit = sqrt(boost::numeric::ublas::norm_2(DataError));
    jiba::rmat InvarSens(ndata, nmod);
    jiba::rvec StartingDataVector(nmeas), DataDiffVector(nmeas);

    std::ofstream misfitfile("misfit.out");
    while (iter < maxiter && gradnorm > mingradnorm && misfit > minmisfit)
      {
        std::cout << "\n\n";
        std::cout << " Iteration: " << iter << std::endl;
        //calculate the response of the starting model
        jiba::ThreeDGravityModel::tTensorMeasVec StartingData(
            Model.CalcTensorGravity());

        for (size_t i = 0; i < ndata; ++i)
          {
            StartingDataVector( i) = CalcInvariant(StartingData.at(i));
            for (size_t j = 0; j < nmod; ++j)
              {
                InvarSens(i, j) = Model.GetTensorSensitivities()(i * 9, j)
                    * StartingData.at(i)(1, 1)
                    + Model.GetTensorSensitivities()(i * 9 + 4, j)
                        * StartingData.at(i)(0, 0)

                + Model.GetTensorSensitivities()(i * 9 + 4, j)
                    * StartingData.at(i)(2, 2)
                    + Model.GetTensorSensitivities()(i * 9 + 8, j)
                        * StartingData.at(i)(1, 1)

                + Model.GetTensorSensitivities()(i * 9, j)
                    * StartingData.at(i)(2, 2)
                    + Model.GetTensorSensitivities()(i * 9 + 8, j)
                        * StartingData.at(i)(0, 0)

                - 2 * Model.GetTensorSensitivities()(i * 9 + 1, j)
                    * StartingData.at(i)(0, 1) - 2
                    * Model.GetTensorSensitivities()(i * 9 + 7, j)
                    * StartingData.at(i)(1, 2) - 2
                    * Model.GetTensorSensitivities()(i * 9 + 2, j)
                    * StartingData.at(i)(0, 2);
              }
          }

        std::transform(DataVector.begin(), DataVector.end(),
            StartingDataVector.begin(), DataDiffVector.begin(), std::minus<
                double>());
        misfit = sqrt(boost::numeric::ublas::norm_2(DataDiffVector));
        std::cout << "Misfit: " << misfit << std::endl;
        misfitfile << iter << " " << misfit << std::endl;
        //depth weighting
        jiba::rvec SensProfile;
        jiba::ExtractMiddleSens(Model, ublas::matrix_range<jiba::rmat>(
            InvarSens, ublas::range(0, nmeas), ublas::range(0, ngrid)), 1,
            SensProfile);

        const double decayexponent = -3.0;
        double z0 = FitZ0(SensProfile, Model.GetZCellSizes(),
            jiba::WeightingTerm(decayexponent));
        std::cout << "Estimated z0: " << z0 << std::endl;

        jiba::rvec WeightVector(zsize), ModelWeight(InvarSens.size2());
        jiba::ConstructDepthWeighting(Model.GetZCellSizes(), z0, WeightVector,
            jiba::WeightingTerm(decayexponent));

        std::fill_n(ModelWeight.begin(), ModelWeight.size(), 0.0);
        for (size_t i = 0; i < nmod; ++i)
          {
            ModelWeight( i) = WeightVector(i % zsize);
          }

        jiba::DataSpaceInversion()(InvarSens, DataDiffVector, ModelWeight,
            DataError, lambda, InvModel);

        gradnorm = sqrt(boost::numeric::ublas::norm_2(InvModel));
        std::cout << "Stepsize: " << gradnorm << std::endl;
        std::cout << "Data Error: " << minmisfit << std::endl;
        //add the result of the inversion to the starting model
        //we only add the gridded part, the  background is always 0 due to the weighting
        std::transform(InvModel.begin(), InvModel.begin() + ngrid,
            Model.SetDensities().origin(), Model.SetDensities().origin(),
            std::plus<double>());
        std::fill(InvModel.begin(), InvModel.end(), 0.0);
        iter++;
      }
    //we want to save the resulting predicted data
    jiba::ThreeDGravityModel::tTensorMeasVec FTGData(Model.CalcTensorGravity());
    for (size_t i = 0; i < ndata; ++i)
      {
        StartingDataVector( i) = CalcInvariant(FTGData.at(i));
      }
    std::transform(DataVector.begin(), DataVector.end(),
        StartingDataVector.begin(), DataDiffVector.begin(),
        std::minus<double>());
    std::transform(DataDiffVector.begin(), DataDiffVector.end(),
        DataError.begin(), DataDiffVector.begin(), std::divides<double>());
    Model.SaveTensorMeasurements(modelfilename + ".inv_ftgdata.nc");
    Model.PlotTensorMeasurements(modelfilename + ".inv_ftgdata.plot");
    Model.WriteVTK(modelfilename + ".inv_ftg.vtk");
    Model.WriteNetCDF(modelfilename + ".inv_ftg.nc");

    std::vector<double> relerror;
    std::copy(DataDiffVector.begin(), DataDiffVector.end(), std::back_inserter(
        relerror));
    jiba::Write3DDataToVTK(modelfilename + ".inv_ftgmisfit.vtk",
        "MisfitInvariant", relerror, Model.GetMeasPosX(), Model.GetMeasPosY(),
        Model.GetMeasPosZ());

    jiba::Write3DTensorDataToVTK(modelfilename + ".inv_ftgdata.vtk", "U",
        FTGData, Model.GetMeasPosX(), Model.GetMeasPosY(), Model.GetMeasPosZ());
  }
