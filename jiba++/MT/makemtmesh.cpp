//============================================================================
// Name        : makegravmesh.cpp
// Author      : Max Moorkamp
// Version     :
// Copyright   : 2008, MM
//============================================================================

/*! \file makemtmesh.cpp
 * Make a netcdf conductivity model file with a specified mesh. The conductivities in this file are all identical to the specified value.
 */

#include <iostream>
#include <string>
#include "../Global/FileUtil.h"
#include "X3DModel.h"

using namespace std;

int main()
  {

    jiba::X3DModel Model;
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

    Model.SetHorizontalCellSize(deltax, deltay, nx, ny);
    Model.SetZCellSizes().resize(boost::extents[nz]);
    Model.SetConductivities().resize(boost::extents[nx][ny][nz]);
    double defaultconductivity = 1.0;
    std::cout << "Conductivity: ";
    std::cin >> defaultconductivity;
    fill_n(Model.SetConductivities().origin(),
        Model.GetConductivities().num_elements(), defaultconductivity);
    fill_n(Model.SetZCellSizes().begin(), nz, deltaz);
    std::string MeshFilename = jiba::AskFilename("Meshfile name: ", false);
    Model.SetFrequencies().assign(1, 10.0);
    std::vector<double> bg_thicknesses(Model.GetZCellSizes().size(), deltaz);
    std::vector<double> bg_conductivities(Model.GetZCellSizes().size(),
        defaultconductivity);
    Model.SetBackgroundThicknesses(bg_thicknesses);
    Model.SetBackgroundConductivities(bg_conductivities);
    Model.WriteNetCDF(MeshFilename);

  }
