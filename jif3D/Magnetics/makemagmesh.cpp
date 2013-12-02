//============================================================================
// Name        : makegravmesh.cpp
// Author      : Max Moorkamp
// Version     :
// Copyright   : 2008, MM
//============================================================================

/*! \file makemagmesh.cpp
 * Make a netcdf magnetic model file with a specified mesh.
 * The susceptibilities in this file are all identical to the specified value.
 */

#include <iostream>
#include <string>
#include "../Global/FileUtil.h"
#include "ThreeDMagneticModel.h"

using namespace std;

int main()
  {

    jif3D::ThreeDMagneticModel Model;
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
    Model.SetMeshSize(nx, ny, nz);

    double defaultsusceptibility = 0.0;
    std::cout << "Susceptibility: ";
    std::cin >> defaultsusceptibility;
    fill_n(Model.SetSusceptibilities().origin(),
        Model.GetSusceptibilities().num_elements(), defaultsusceptibility);
    fill_n(Model.SetXCellSizes().begin(), nx, deltax);
    fill_n(Model.SetYCellSizes().begin(), ny, deltay);
    fill_n(Model.SetZCellSizes().begin(), nz, deltaz);
    std::string MeshFilename = jif3D::AskFilename("Meshfile name: ", false);
    Model.WriteNetCDF(MeshFilename);

  }
