//============================================================================
// Name        : gravinv.cpp
// Author      : Max Moorkamp
// Version     :
// Copyright   : 2008, MM
//============================================================================


#include <iostream>
#include <fstream>
#include <string>
#include "ThreeDGravityModel.h"
#include "../Inversion/MatrixTools.h"
#include "ReadWriteGravityData.h"
#include <boost/numeric/bindings/atlas/cblas2.hpp>
#include <boost/numeric/bindings/atlas/cblas3.hpp>

namespace atlas = boost::numeric::bindings::atlas;
namespace ublas = boost::numeric::ublas;
int main(int argc, char *argv[])
  {
    jiba::ThreeDGravityModel Model(true);

    jiba::ThreeDGravityModel::tScalarMeasVec Data;
    jiba::ThreeDGravityModel::tMeasPosVec PosX, PosY, PosZ;

    std::string modelfilename, datafilename;
    std::cout << "Mesh Filename: ";
    std::cin >> modelfilename;
    //we read in a complete modelfile, but we only use the mesh information
    Model.ReadNetCDF(modelfilename);
    //get the name of the file containing the data and read it in
    std::cout << "Data Filename: ";
    std::cin >> datafilename;
    jiba::ReadScalarGravityMeasurements(datafilename, Data, PosX, PosY, PosZ);
    const size_t nmeas = PosX.size();
    //set the measurement points in the model to those of the data
    Model.ClearMeasurementPoints();
    for (size_t i = 0; i < nmeas; ++i)
      {
        Model.AddMeasurementPoint(PosX.at(i), PosY.at(i), PosZ.at(i));
      }
    //calculate the response, we don't actually care about the densities
    //and the model response, but we need the sensitivity matrix
    Model.CalcGravity();
    Model.SaveScalarMeasurements(modelfilename + ".new.nc");
    jiba::rmat Sensitivities(Model.GetScalarSensitivities());

    double z0;
    std::cout << "Enter z0: ";
    std::cin >> z0;

    const size_t xsize = Model.GetDensities().shape()[0];
    const size_t ysize = Model.GetDensities().shape()[1];
    const size_t zsize = Model.GetDensities().shape()[2];
    const size_t nmod = xsize * ysize * zsize;
    jiba::rvec WeightVector(zsize);
    jiba::ConstructDepthWeighting(Model.GetXCellSizes(), Model.GetYCellSizes(),
        Model.GetZCellSizes(), z0, WeightVector);

    jiba::rvec ModelCovar(nmod);
    ModelCovar *= 0.0;
    std::copy(WeightVector.begin(), WeightVector.end(), ModelCovar.begin());
    jiba::rmat FilteredSens(nmod, nmeas);
    for (size_t i = 0; i < nmod; ++i)
      {
        boost::numeric::ublas::matrix_row<jiba::rmat>
            CurrentRow(FilteredSens, i);

        atlas::gemv(CblasNoTrans, 1.0, Sensitivities, ModelCovar, 0.0,
            CurrentRow);
        std::rotate(ModelCovar.begin(), ModelCovar.end() - zsize,
            ModelCovar.end());
      }

    jiba::rmat Gamma(nmeas, nmeas);
    atlas::gemm(CblasNoTrans, CblasNoTrans, 1.0, Sensitivities, FilteredSens,
        0.0, Gamma);
    //we can play around with the threshold for the included eigenvalues
    double upthresh, lowthresh;
    std::cout << "Upper threshold: ";
    std::cin >> upthresh;
    std::cout << "Lower threshold: ";
    std::cin >> lowthresh;

    jiba::rmat DataInverse(nmeas, nmeas);
    jiba::GeneralizedInverse(Gamma, DataInverse, lowthresh, upthresh);

    jiba::rvec InvModel(xsize * ysize * zsize);

    jiba::rmat ModelInverse(nmod, nmeas);
    atlas::gemm(CblasNoTrans, CblasNoTrans, 1.0, FilteredSens, DataInverse,
        0.0, ModelInverse);

    atlas::gemv(ModelInverse, Data, InvModel);

    std::copy(InvModel.begin(), InvModel.end(), Model.SetDensities().origin());
    Model.CalcGravity();
    Model.SaveScalarMeasurements(modelfilename + ".inv_data.nc");
    Model.PlotScalarMeasurements(modelfilename + ".inv.plot");
    Model.WriteVTK(modelfilename + ".inv.vtk");
    Model.WriteNetCDF(modelfilename + ".inv.nc");
  }
