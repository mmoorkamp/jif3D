//============================================================================
// Name        : makedcresmesh.cpp
// Author      : Zhanjie Shi and Richard.W Hobbs
// Version     : April 2014
// Copyright   : 2014, Zhanjie Shi and Richard.W Hobbs
//============================================================================

/*! \file makedcresmesh.cpp
 * Make a netcdf resistivity model file with a specified mesh. The resistivities in this file are all identical to the specified value.
 */

#include <iostream>
#include <string>
#include "../Global/FileUtil.h"
#include "ThreeDDCResistivityModel.h"

using namespace std;

/*! \file makedcresmesh.cpp
 * Make a netcdf resistivity model file with a specified mesh. The program asks for the size in the three coordinate directions and
 * the cell size for each direction. All cells will have the same size for x and y direction and have varying cell sizes in z-direction.
 * The program also asks for a resistivity value that the mesh will be filled with.
 */

int main()
  {

    jif3D::ThreeDDCResistivityModel Model;
    int nx, ny, nz;
    double deltax, deltay;
    vector<double> deltaz;
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
    deltaz.resize(nz);
    for (int i = 0; i < nz; i++)
      {
        cin >> deltaz[i];
      }

    //set the cell sizes and allocate memory for the mesh
    Model.SetHorizontalCellSize(deltax, deltay, nx, ny);
    Model.SetZCellSizes().resize(boost::extents[nz]);
    for (int j = 0; j < nz; j++)
      {
        Model.SetZCellSizes()[j] = deltaz[j];
      }
    Model.SetResistivities().resize(boost::extents[nx][ny][nz]);
    //ask for a resistivity to fill the mesh with
    double defaultresistivity = 1.0;
    std::cout << "Resistivity: ";
    std::cin >> defaultresistivity;
    fill_n(Model.SetResistivities().origin(), Model.GetResistivities().num_elements(),
        defaultresistivity);
    //ask for a filename to write the mesh to
    std::string MeshFilename = jif3D::AskFilename("Meshfile name: ", false);
    Model.WriteNetCDF(MeshFilename);

  }
