//============================================================================
// Name        : makegravmesh.cpp
// Author      : Max Moorkamp
// Version     :
// Copyright   : 2008, MM
//============================================================================

/*! \file makegravmesh.cpp
 * Make a netcdf gravity model file with a specified mesh.
 * The densities in this file are all identical to the specified value.
 */

#include "ThreeDGravityModel.h"
#include "../Global/FileUtil.h"
#include <iostream>
#include <string>

using namespace std;

int main()
  {

    jif3D::ThreeDGravityModel Model;
    int nx, ny, nz;
    double deltax, deltay, deltaz;
    cout << "Number of cells in x-direction: ";
    cin >> nx;
    cout << "Number of cells in y-direction: ";
    cin >> ny;
    cout << "Number of cells in z-direction: ";
    cin >> nz;
    cout << "Cells size in x-direction (m): ";
    cin >> deltax;
    cout << "Cells size in y-direction (m): ";
    cin >> deltay;
    cout << "Cells size in z-direction (m): ";
    cin >> deltaz;

    Model.SetMeshSize(nx, ny, nz);
    double defaultdensity = 0.0;
    std::cout << "Density: ";
    std::cin >> defaultdensity;
    fill_n(Model.SetDensities().origin(), Model.GetDensities().num_elements(),
        defaultdensity);
    jif3D::ThreeDModelBase::t3DModelDim XCS(nx, deltax), YCS(ny, deltay), ZCS(nz, deltaz);
    Model.SetXCellSizes(XCS);
    Model.SetYCellSizes(YCS);
    Model.SetZCellSizes(ZCS);

    std::vector<double> bg_densities(nz, defaultdensity), bg_thicknesses(nz, deltaz);
    Model.SetBackgroundDensities(bg_densities);
    Model.SetBackgroundThicknesses(bg_thicknesses);
    std::string MeshFilename = jif3D::AskFilename("Meshfile name: ", false);
    Model.WriteNetCDF(MeshFilename);
    Model.WriteVTK(MeshFilename + ".vtk");
  }
