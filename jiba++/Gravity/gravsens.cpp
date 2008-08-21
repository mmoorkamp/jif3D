//============================================================================
// Name        : gravsens.cpp
// Author      : Max Moorkamp
// Version     :
// Copyright   : 2008, MM
//============================================================================


#include <iostream>
#include <fstream>
#include <string>
#include "ThreeDGravityModel.h"
#include "ReadWriteGravityData.h"
#include "../Inversion/MatrixTools.h"
#include <boost/numeric/bindings/atlas/cblas2.hpp>
#include "../ModelBase/VTKTools.h"
namespace atlas = boost::numeric::bindings::atlas;

int main(int argc, char *argv[])
  {
    jiba::ThreeDGravityModel Model(true);

    std::string modelfilename, datafilename;
    std::cout << "Mesh Filename: ";
    std::cin >> modelfilename;
    //we read in a complete modelfile, but we only use the mesh information
    Model.ReadNetCDF(modelfilename);

    //we also read in data, but we only use the information about the measurements
    jiba::ThreeDGravityModel::tScalarMeasVec Data;
    jiba::ThreeDGravityModel::tMeasPosVec PosX, PosY, PosZ;
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
    jiba::rmat Sensitivities(Model.GetScalarSensitivities());

    jiba::rvec s;
    jiba::rmat u, vt;
    jiba::SVD(Sensitivities, s, u, vt);
    std::ofstream evalfile((modelfilename + ".eval").c_str());
    std::copy(s.begin(), s.end(), std::ostream_iterator<double>(evalfile, "\n"));
    std::cout << "There are " << s.size() << " eigenvalues.\n";
    int index = 0;
    std::cout << "Enter an index for the eigenvalue to plot: ";
    std::cin >> index;
    const size_t xsize = Model.GetDensities().shape()[0];
    const size_t ysize = Model.GetDensities().shape()[1];
    const size_t zsize = Model.GetDensities().shape()[2];
    jiba::ThreeDModelBase::t3DModelData sens(boost::extents[xsize][ysize][zsize]);
    std::copy(row(vt, index).begin(),row(vt, index).end(),sens.data());
    jiba::Write3DModelToVTK(modelfilename + ".sens.vtk", "scalar_sens", Model.GetXCellSizes(),
        Model.GetYCellSizes(), Model.GetZCellSizes(), sens);
  }
