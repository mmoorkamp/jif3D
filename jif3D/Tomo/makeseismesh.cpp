//============================================================================
// Name        : makeseismesh.cpp
// Author      : Max Moorkamp
// Version     :
// Copyright   : 2009, MM
//============================================================================

/*! \file makeseismesh.cpp
 * Make a netcdf seismic model file with a specified mesh. We can specify the velocity in the top layer
 * and in the bottom layer in m/s. The velocities in the other grid cells will be linearly
 * interpolated with depth.
 */

#include <iostream>
#include <string>

#include "../Tomo/ThreeDSeismicModel.h"
#include "../Global/FileUtil.h"
using namespace std;

int main()
  {
    jif3D::ThreeDSeismicModel Model;
    int nx, ny, nz;
    double cellsize;
    //first we enter the number of cells in each direction
    cout << "Number of cells in x-direction: ";
    cin >> nx;
    cout << "Number of cells in y-direction: ";
    cin >> ny;
    cout << "Number of cells in z-direction: ";
    cin >> nz;
    //for the forward code all cells have to have the same size in each direction
    //so we only need to ask for one cell size
    cout << "Cell size in m: ";
    cin >> cellsize;

    Model.SetCellSize(cellsize, nx, ny, nz);
    //we need a velocity gradient for sensible seismic refraction calculation
    //so we ask for the velocity at the top of the model
    //and at the bottom of the model and fill each layer with a linear interpolation
    //between the two values.
    double topvel = 0.0;
    std::cout << "Velocity at top in m/s: ";
    std::cin >> topvel;
    double bottomvel = 0.0;
    std::cout << "Velocity at bottom  in m/s: ";
    std::cin >> bottomvel;
    const double firstdepth = Model.GetZCoordinates()[0];
    const double bottomdepth = Model.GetZCoordinates()[nz - 1];
    //calculate the depth for each cell and the corresponding velocity
    //then assign the resulting slowness to the cell
    for (size_t i = 0; i < Model.GetSlownesses().num_elements(); ++i)
      {
        double Depth = Model.GetZCoordinates()[i % nz];
        double Velocity = topvel
            + (Depth - firstdepth) * (bottomvel - topvel) / (bottomdepth - firstdepth);
        Model.SetSlownesses().origin()[i] = 1.0 / Velocity;
      }
    //finally we ask for a name for the meshfile to write to
    //and save the created mesh with the velocity information
    std::string MeshFilename = jif3D::AskFilename("Meshfile name: ", false);
    Model.WriteNetCDF(MeshFilename);
    Model.WriteVTK(MeshFilename + ".vtk");
  }
