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

namespace atlas = boost::numeric::bindings::atlas;

int main(int argc, char *argv[])
  {
    jiba::ThreeDGravityModel Model(true);

    jiba::ThreeDGravityModel::tScalarMeasVec Data;
    jiba::ThreeDGravityModel::tMeasPosVec PosX, PosY, PosZ;

    std::string modelfilename, datafilename;
    std::cout << "Model Filename: ";
    std::cin >> modelfilename;
    Model.ReadNetCDF(modelfilename);
    std::cout << "Data Filename: ";
    std::cin >> datafilename;
    jiba::ReadScalarGravityMeasurements(datafilename, Data, PosX, PosY, PosZ);
    const size_t nmeas = PosX.size();
    Model.ClearMeasurementPoints();
    for (size_t i = 0; i < nmeas; ++i)
      {
        Model.AddMeasurementPoint(PosX.at(i), PosY.at(i), PosZ.at(i));
      }
    Model.CalcGravity();
    Model.SaveScalarMeasurements(modelfilename + ".new.nc");
    jiba::rmat Sensitivities(Model.GetScalarSensitivities());

    double upthresh, lowthresh;
    std::cout << "Upper threshold: ";
    std::cin >> upthresh;
    std::cout << "Lower threshold: ";
    std::cin >> lowthresh;
    jiba::rmat Inverse(Sensitivities.size2(), Sensitivities.size1());
    jiba::GeneralizedInverse(Sensitivities, Inverse, lowthresh, upthresh);
    const size_t xsize = Model.GetDensities().shape()[0];
    const size_t ysize = Model.GetDensities().shape()[1];
    const size_t zsize = Model.GetDensities().shape()[2];
    jiba::rvec InvModel(xsize * ysize * zsize);
    atlas::gemv(Inverse, Data, InvModel);

    for (size_t i = 0; i < xsize; ++i)
      {
        for (size_t j = 0; j < ysize; ++j)
          {
            for (size_t k = 0; k < zsize; ++k)
              {
                Model.SetDensities()[i][j][k] = InvModel(i * (ysize * zsize)
                    + j * zsize + k);
              }
          }
      }
    Model.SaveScalarMeasurements(modelfilename + ".inv_data.nc");
    Model.PlotScalarMeasurements(modelfilename + ".inv.plot");
    Model.WriteVTK(modelfilename + ".inv.vtk");
    Model.WriteNetCDF(modelfilename + ".inv.nc");
  }
