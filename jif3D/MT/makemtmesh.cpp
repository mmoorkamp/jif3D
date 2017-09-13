//============================================================================
// Name        : makegravmesh.cpp
// Author      : Max Moorkamp
// Version     :
// Copyright   : 2008, MM
//============================================================================

/*! \file makemtmesh.cpp
 * Make a netcdf conductivity model file with a specified mesh. The conductivities in this file are all identical to the specified value.
 */

#include "../Global/Jif3DGlobal.h"
#include "../Global/FileUtil.h"
#include "X3DModel.h"
#include <boost/program_options.hpp>
#include <iostream>
#include <string>


using namespace std;

/*! \file makemtmesh.cpp
 * Make a forward modeling mesh for X3D. The program asks for the size in the three coordinate directions and
 * the cell size for each direction. All cells will have the same size for the two horizontal directions.
 * For the vertical direction we can specify the thickness of the top layer and a factor by which each cell
 * increases in thickness compared to the previous layer.
 * The program also asks for a conductivity value that the mesh will be filled with.
 * It will also create a layered background with the same number of layers as cells in x-direction and filled with
 * the same conductivity value as the mesh.
 */
namespace po = boost::program_options;

int main(int argc, char *argv[])
  {
    int incstart = 0;
    double rounding = 1.0;
    po::options_description desc("General options");
    desc.add_options()("help", "produce help message")("rounding",
        po::value<double>(&rounding),
        "Round layer thicknesses to multiple of this number in meters.")("incstart",
        po::value<int>(&incstart),
        "Index of the layer where to start increasing the cell size in z-direction");

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);

    jif3D::X3DModel Model;
    int nx, ny, nz;
    double deltax, deltay, deltaz;
    //first find out the basic mesh parameters
    //the number of cells in each coordinate direction
    cout << "Number of cells in x-direction:: ";
    cin >> nx;
    cout << "Number of cells in y-direction:: ";
    cin >> ny;
    cout << "Number of cells in z-direction:: ";
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
    Model.SetMeshSize(nx, ny, nz);
    Model.SetHorizontalCellSize(deltax, deltay, nx, ny);
    //the size of each cell in z-direction increases by the specified factor for each layer
    for (int i = 0; i < nz; ++i)
      {
        double thickness = deltaz;
        if (i >= incstart)
          thickness *= pow(factor, i-incstart);
        //x3d has some problems handling thicknesses over 10km with full meter precision
        //so if the thickness is > 10km we round to 100m
        thickness = floor(thickness / rounding) * rounding;
        Model.SetZCellSizes()[i] = thickness;
      }

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
    std::vector<double> bg_thicknesses(Model.GetZCellSizes().size());
    std::copy(Model.GetZCellSizes().begin(), Model.GetZCellSizes().end(),
        bg_thicknesses.begin());
    std::vector<double> bg_conductivities(Model.GetZCellSizes().size(),
        defaultconductivity);
    Model.SetBackgroundThicknesses(bg_thicknesses);
    Model.SetBackgroundConductivities(bg_conductivities);
    Model.WriteNetCDF(MeshFilename);
    Model.WriteVTK(MeshFilename + ".vtk");
    Model.WriteModEM(MeshFilename + ".dat");

  }
