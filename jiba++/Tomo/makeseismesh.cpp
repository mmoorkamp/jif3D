//============================================================================
// Name        : makeseismesh.cpp
// Author      : Max Moorkamp
// Version     :
// Copyright   : 2009, MM
//============================================================================

/*! \file makeseismesh.cpp
 * Make a netcdf seismic model file with a specified mesh. We can specify the velocity in the top layer
 * and in the bottom lauer in m/s. The velocities in the other grid cells will be linearly
 * interpolated with depth.
 */

#include <iostream>
#include <string>
#include "ThreeDSeismicModel.h"
#include "../Global/FileUtil.h"
using namespace std;

int main(int argc, char *argv[])
  {
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

    Model.SetCellSize(cellsize, nx, ny, nz);
    Model.SetSlownesses().resize(boost::extents[nx][ny][nz]);
    double topvel = 0.0;
    std::cout << "Velocity at top: ";
    std::cin >> topvel;
    double bottomvel = 0.0;
    std::cout << "Velocity at bottom: ";
    std::cin >> bottomvel;
    const double firstdepth = Model.GetZCoordinates()[0];
    const double bottomdepth = Model.GetZCoordinates()[nz -1];
    for (size_t i = 0; i < Model.GetSlownesses().num_elements(); ++i)
      {
        double Depth = Model.GetZCoordinates()[i % nz];
        double Velocity = topvel + (Depth - firstdepth) * (bottomvel - topvel)/(bottomdepth - firstdepth);
        Model.SetSlownesses().origin()[i] = 1.0/Velocity;
      }


    std::string MeshFilename = jiba::AskFilename("Meshfile name: ", false);
    Model.WriteNetCDF(MeshFilename);
  }
