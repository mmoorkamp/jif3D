//============================================================================
// Name        : ftginv.cpp
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
    jiba::rvec DataVector(Data.size() * 9);
    for (size_t i = 0; i < Data.size(); ++i)
      {
        copy(Data.at(i).data().begin(),
            Data.at(i).data().end(), DataVector.begin()+i*9);
      }

    const size_t nmeas = PosX.size();
    //set the measurement points in the model to those of the data
    Model.ClearMeasurementPoints();
    for (size_t i = 0; i < nmeas; ++i)
      {
        Model.AddMeasurementPoint(PosX.at(i), PosY.at(i), PosZ.at(i));
      }
    //calculate the response, we don't actually care about the densities
    //and the model response, but we need the sensitivity matrix
    Model.CalcTensorGravity();
    Model.SaveTensorMeasurements(modelfilename + ".new.nc");
    jiba::rmat Sensitivities(Model.GetTensorSensitivities());

    //we can play around with the threshold for the included eigenvalues
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
    atlas::gemv(Inverse, DataVector, InvModel);

    std::copy(InvModel.begin(), InvModel.end(), Model.SetDensities().origin());
    Model.CalcTensorGravity();
    Model.SaveTensorMeasurements(modelfilename + ".inv_data.nc");
    Model.PlotTensorMeasurements(modelfilename + ".inv.plot");
    Model.WriteVTK(modelfilename + ".inv.vtk");
    Model.WriteNetCDF(modelfilename + ".inv.nc");
  }
