//============================================================================
// Name        : makeseismesh.cpp
// Author      : Max Moorkamp
// Version     :
// Copyright   : 2009, MM
//============================================================================

/*! \file makeseismesh.cpp
 * Make a netcdf seismic model file with a specified mesh. The slownesses in this file are all identical to the specified value.
 */


#include <iostream>
#include <string>
#include "ThreeDSeismicModel.h"

using namespace std;

int main(int argc, char *argv[])
  {
    std::string MeshFilename;
    jiba::ThreeDSeismicModel Model;
    int nx, ny, nz;
    double cellsize;
    cout << "Nx: ";
    cin >> nx;
    cout << "Ny: ";
    cin >> ny;
    cout << "Nz: ";
    cin >> nz;
    cout << "Cell size: ";
    cin >> cellsize;


    Model.SetCellSize(cellsize,nz,ny,nz);
    Model.SetSlownesses().resize(boost::extents[nx][ny][nz]);
    double defaultslowness = 0.0;
    std::cout << "Slowness: ";
    std::cin >> defaultslowness;
    fill_n(Model.SetSlownesses().origin(),Model.GetSlownesses().num_elements(),defaultslowness);
    cout << "Meshfile name: ";
    cin >> MeshFilename;
    Model.WriteNetCDF(MeshFilename);


  }
