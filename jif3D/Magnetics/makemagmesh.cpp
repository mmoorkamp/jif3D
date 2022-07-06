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

#include "../Global/Serialization.h"
#include "../Global/FileUtil.h"
#include "ThreeDSusceptibilityModel.h"
#include <iostream>
#include <string>

using namespace std;

int main()
  {

    jif3D::ThreeDSusceptibilityModel Model;
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
    jif3D::ThreeDModelBase::t3DModelDim XCS(nx, deltax), YCS(ny, deltay), ZCS(nz, deltaz);
    Model.SetXCellSizes(XCS);
    Model.SetYCellSizes(YCS);
    Model.SetZCellSizes(ZCS);
    std::string MeshFilename = jif3D::AskFilename("Meshfile name: ", false);
    Model.WriteNetCDF(MeshFilename);

  }
