//============================================================================
// Name        : makegravmesh.cpp
// Author      : Max Moorkamp
// Version     :
// Copyright   : 2008, MM
//============================================================================

/*! \file makegravmesh.cpp
 * Make a netcdf gravity model file with a specified mesh. The densities in this file are all identical to the specified value.
 */


#include <iostream>
#include <string>
#include "../Global/FileUtil.h"
#include "ThreeDGravityModel.h"

using namespace std;

int main()
  {

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
    double defaultdensity = 0.0;
    std::cout << "Density: ";
    std::cin >> defaultdensity;
    fill_n(Model.SetDensities().origin(),Model.GetDensities().num_elements(),defaultdensity);
    fill_n(Model.SetXCellSizes().begin(),nx,deltax);
    fill_n(Model.SetYCellSizes().begin(),ny,deltay);
    fill_n(Model.SetZCellSizes().begin(),nz,deltaz);
    std::string MeshFilename = jiba::AskFilename("Meshfile name: ",false);
    Model.WriteNetCDF(MeshFilename);


  }
