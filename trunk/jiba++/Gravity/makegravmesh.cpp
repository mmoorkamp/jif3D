//============================================================================
// Name        : makegravmesh.cpp
// Author      : Max Moorkamp
// Version     :
// Copyright   : 2008, MM
//============================================================================


#include <iostream>
#include <string>
#include "ThreeDGravityModel.h"

using namespace std;

int main(int argc, char *argv[])
  {
    std::string MeshFilename;
    jiba::ThreeDGravityModel Model;
    int nx, ny, nz;
    double deltax, deltay, deltaz;
    cout << "Nx: ";
    cin >> nx;
    cout << "Ny: ";
    cin >> ny;
    cout << "Nz: ";
    cin >> nz;
    cout << "Delta x: ";
    cin >> deltax;
    cout << "Delta y: ";
    cin >> deltay;
    cout << "Delta z: ";
    cin >> deltaz;

    Model.SetXCellSizes().resize(boost::extents[nx]);
    Model.SetYCellSizes().resize(boost::extents[ny]);
    Model.SetZCellSizes().resize(boost::extents[nz]);
    Model.SetDensities().resize(boost::extents[nx][ny][nz]);

    fill_n(Model.SetXCellSizes().begin(),nx,deltax);
    fill_n(Model.SetYCellSizes().begin(),ny,deltay);
    fill_n(Model.SetZCellSizes().begin(),nz,deltaz);
    cout << "Meshfile name: ";
    cin >> MeshFilename;
    Model.WriteNetCDF(MeshFilename);


  }
