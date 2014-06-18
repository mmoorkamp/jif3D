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

/*! \file makemtmesh.cpp
 * Make a forward modeling mesh for X3D. The program asks for the size in the three coordinate directions and
 * the cell size for each direction. All cells will have the same size for one direction even though x3d supports
 * varying cell sizes in z-direction. The program also asks for a conductivity value that the mesh will be filled with.
 * It will also create a layered background with the same number of layers as cells in x-direction and filled with
 * the same conductivity value as the mesh.
 */

int main()
  {

    jif3D::X3DModel Model;
    int nx, ny, nz;
    double deltax, deltay, deltaz;
    //first find out the basic mesh parameters
    //the number of cells in each coordinate direction
    cout << "Nx: ";
    cin >> nx;
    cout << "Ny: ";
    cin >> ny;
    cout << "Nz: ";
    cin >> nz;
    //and the size of the cells in each direction
    cout << "Cell size x [m]: ";
    cin >> deltax;
    cout << "Cell size y [m]: ";
    cin >> deltay;
    cout << "Cell size z [m]: ";
    cin >> deltaz;

    double factor = 1.0;
    cout << "Increase factor for z: ";
    cin >> factor;
    //set the cell sizes and allocate memory for the mesh
    Model.SetHorizontalCellSize(deltax, deltay, nx, ny);
    Model.SetZCellSizes().resize(boost::extents[nz]);
    for (size_t i = 0; i < nz; ++i)
      {
        Model.SetZCellSizes()[i] = floor(deltaz * pow(factor, i));

      }
    Model.SetMeshSize(nx, ny, nz);

    //ask for a conductivity to fill the mesh with
    double defaultconductivity = 1.0;
    std::cout << "Conductivity: ";
    std::cin >> defaultconductivity;
    fill_n(Model.SetConductivities().origin(), Model.GetConductivities().num_elements(),
        defaultconductivity);
    //ask for a filename to write the mesh to
    std::string MeshFilename = jif3D::AskFilename("Meshfile name: ", false);
    //we set the calculation frequencies to a dummy value
    //the forward modeling program asks for those anyway
    Model.SetFrequencies().assign(1, 10.0);
    //fill the background
    std::vector<double> bg_thicknesses(Model.GetZCellSizes().size(), deltaz);
    std::vector<double> bg_conductivities(Model.GetZCellSizes().size(),
        defaultconductivity);
    Model.SetBackgroundThicknesses(bg_thicknesses);
    Model.SetBackgroundConductivities(bg_conductivities);
    Model.WriteNetCDF(MeshFilename);

  }
