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
#include <cstdlib>
#include <omp.h>
#include <boost/program_options.hpp>
#include "../Global/FileUtil.h"
#include "../Global/convert.h"
#include "../MT/X3DModel.h"
#include "../MT/X3DMTCalculator.h"
#include "../MT/ReadWriteImpedances.h"

using namespace std;
namespace po = boost::program_options;
/*! \file makemtmesh.cpp
 * Make a forward modeling mesh for X3D. The program asks for the size in the three coordinate directions and
 * the cell size for each direction. All cells will have the same size for one direction even though x3d supports
 * varying cell sizes in z-direction. The program also asks for a conductivity value that the mesh will be filled with.
 * It will also create a layered background with the same number of layers as cells in x-direction and filled with
 * the same conductivity value as the mesh.
 */

int main(int argc, char *argv[])
  {

    po::options_description desc("General options");
    desc.add_options()("help", "produce help message")("threads", po::value<int>(),
        "The number of openmp threads")("topthick",
        "Thickness of the top layer in m");
    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);

    if (vm.count("help"))
      {
        std::cout << desc << "\n";
        return 1;
      }
    if (vm.count("threads"))
      {
        omp_set_num_threads(vm["threads"].as<int>());
      }

    jif3D::X3DModel Model;
    size_t nx, ny, nz;
    double deltax, deltay, deltaz;
    //first find out the basic mesh parameters
    //the number of cells in each coordinate direction
    cout << "Number of cells Nx: ";
    cin >> nx;
    cout << "Number of cells Ny: ";
    cin >> ny;
    cout << "Number of cells Nz: ";
    cin >> nz;
    //and the size of the cells in each direction
    cout << "Cell size x [m]: ";
    cin >> deltax;
    cout << "Cell size y [m]: ";
    cin >> deltay;
    cout << "Cell size z [m]: ";
    cin >> deltaz;
    //set the cell sizes and allocate memory for the mesh
    Model.SetHorizontalCellSize(deltax, deltay, nx, ny);
    Model.SetZCellSizes().resize(boost::extents[nz]);
    fill_n(Model.SetZCellSizes().begin(), nz, deltaz);
    if (vm.count("topthick"))
      {
        Model.SetZCellSizes()[0] = vm["topthick"].as<double>();
      }
    Model.SetConductivities().resize(boost::extents[nx][ny][nz]);
    //ask for a conductivity to fill the mesh with
    double bg_conductivity = 1.0;
    std::cout << "Background Conductivity [S/m] : ";
    std::cin >> bg_conductivity;
    double phase1cond, phase2cond, phase1frac;
    std::cout << "Conductivity Phase 1 [S/m] : ";
    std::cin >> phase1cond;
    std::cout << "Conductivity Phase 2 [S/m] : ";
    std::cin >> phase2cond;
    std::cout << "Fraction Phase 1: ";
    std::cin >> phase1frac;
    srand48(time(0));
    double frequency;
    std::cout << "Frequency [Hz] : ";
    std::cin >> frequency;
    double posx, posy, posz = 0;
    posx = (deltax * nx) / 2.0;
    posy = (deltay * ny) / 2.0;
    Model.SetFrequencies().assign(1, frequency);
    Model.AddMeasurementPoint(posx, posy, posz);

    //ask for a filename to write the mesh to
    std::string OutFilename = jif3D::AskFilename("Outfile name: ", false);
    //we set the calculation frequencies to a dummy value
    //the forward modeling program asks for those anyway

    //fill the background
    std::vector<double> bg_thicknesses(Model.GetZCellSizes().size(), deltaz);
    std::vector<double> bg_conductivities(Model.GetZCellSizes().size(), bg_conductivity);
    Model.SetBackgroundThicknesses(bg_thicknesses);
    Model.SetBackgroundConductivities(bg_conductivities);
    size_t nrealmax;
    std::cout << "Realizations: ";
    std::cin >> nrealmax;

    std::ofstream zxy((OutFilename + "_zxy.out").c_str());
    std::ofstream zyx((OutFilename + "_zyx.out").c_str());
    std::ofstream zxx((OutFilename + "_zxx.out").c_str());
    std::ofstream zyy((OutFilename + "_zyy.out").c_str());

#pragma omp parallel for shared(Model, zxy, zyx, zxx, zyy)
    for (size_t nreal = 0; nreal < nrealmax; ++nreal)
      {
        jif3D::X3DModel RealModel(Model);
        for (size_t i = 0; i < nx; ++i)
          {
            for (size_t j = 0; j < ny; ++j)
              {
                RealModel.SetConductivities()[i][j][0] = bg_conductivity;
                for (size_t k = 1; k < nz; ++k)
                  {
                    if (drand48() < phase1frac)
                      {
                        RealModel.SetConductivities()[i][j][k] = phase1cond;
                      }
                    else
                      {
                        RealModel.SetConductivities()[i][j][k] = phase2cond;
                      }
                  }
              }
          }
        std::string realstring(jif3D::stringify(nreal));

        jif3D::X3DMTCalculator Calculator;
        jif3D::rvec Impedances(Calculator.Calculate(RealModel));
        jif3D::rvec Errors(Impedances.size(), 0.0);
#pragma omp critical(write_files)
          {
            std::cout << "Realization: " << nreal << std::endl;
            RealModel.WriteNetCDF(OutFilename + realstring + ".nc");
            RealModel.WriteVTK(OutFilename + realstring + ".vtk");
            zxx << nreal << " " << Impedances(0) << " " << Impedances(1) << std::endl;
            zxy << nreal << " " << Impedances(2) << " " << Impedances(3) << std::endl;
            zyx << nreal << " " << Impedances(4) << " " << Impedances(5) << std::endl;
            zyy << nreal << " " << Impedances(6) << " " << Impedances(7) << std::endl;
          }
      }

  }
