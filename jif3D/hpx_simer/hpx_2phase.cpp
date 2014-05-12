//============================================================================
// Name        : makegravmesh.cpp
// Author      : Max Moorkamp
// Version     :
// Copyright   : 2008, MM
//============================================================================

#include <hpx/config.hpp>
#include <hpx/hpx_init.hpp>
#include <hpx/include/actions.hpp>
#include <hpx/include/util.hpp>
#include <hpx/include/lcos.hpp>
#include <hpx/include/iostreams.hpp>

#include <iostream>
#include <string>
#include <cstdlib>

#include "../Global/FileUtil.h"
#include "../Global/convert.h"
#include "../MT/X3DModel.h"
#include "../MT/X3DMTCalculator.h"

jif3D::rvec CalcRealization(jif3D::X3DModel Model, double bg_conductivity,
    double phase1cond, double phase2cond, double phase1frac);
HPX_PLAIN_ACTION(CalcRealization, Calc_action);

jif3D::rvec CalcRealization(jif3D::X3DModel Model, double bg_conductivity,
    double phase1cond, double phase2cond, double phase1frac)
  {
    const size_t nx = Model.GetXCoordinates().num_elements();
    const size_t ny = Model.GetYCoordinates().num_elements();
    const size_t nz = Model.GetZCoordinates().num_elements();

    for (size_t i = 0; i < nx; ++i)
      {
        for (size_t j = 0; j < ny; ++j)
          {
            Model.SetConductivities()[i][j][0] = bg_conductivity;
            for (size_t k = 1; k < nz; ++k)
              {
                if (drand48() < phase1frac)
                  {
                    Model.SetConductivities()[i][j][k] = phase1cond;
                  }
                else
                  {
                    Model.SetConductivities()[i][j][k] = phase2cond;
                  }
              }
          }
      }

    jif3D::X3DMTCalculator Calculator;
    jif3D::rvec Impedances(Calculator.Calculate(Model));
    return Impedances;
  }

int hpx_main(int argc, char* argv[])
  {

    using hpx::cout;
    using std::cin;

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
    Model.SetMeshSize(nx, ny, nz);
    Model.SetHorizontalCellSize(deltax, deltay, nx, ny);
    std::fill_n(Model.SetZCellSizes().begin(), nz, deltaz);
    //ask for a conductivity to fill the mesh with
    double bg_conductivity = 1.0;
    cout << "Background Conductivity [S/m] : ";
    cin >> bg_conductivity;
    double phase1cond, phase2cond, phase1frac;
    cout << "Conductivity Phase 1 [S/m] : ";
    cin >> phase1cond;
    cout << "Conductivity Phase 2 [S/m] : ";
    cin >> phase2cond;
    cout << "Fraction Phase 1: ";
    cin >> phase1frac;
    srand48(time(0));
    double frequency;
    cout << "Frequency [Hz] : ";
    cin >> frequency;
    double posx, posy, posz = 0;
    posx = (deltax * nx) / 2.0;
    posy = (deltay * ny) / 2.0;
    Model.SetFrequencies().assign(1, frequency);
    Model.AddMeasurementPoint(posx, posy, posz);

    //ask for a filename to write the mesh to
    std::string OutFilename = jif3D::AskFilename("Outfile name: ", false);

    //fill the background
    std::vector<double> bg_thicknesses(Model.GetZCellSizes().size(), deltaz);
    std::vector<double> bg_conductivities(Model.GetZCellSizes().size(), bg_conductivity);
    Model.SetBackgroundThicknesses(bg_thicknesses);
    Model.SetBackgroundConductivities(bg_conductivities);
    size_t nrealmax;
    cout << "Realizations: ";
    cin >> nrealmax;

    std::ofstream zxy((OutFilename + "_zxy.out").c_str());
    std::ofstream zyx((OutFilename + "_zyx.out").c_str());
    std::ofstream zxx((OutFilename + "_zxx.out").c_str());
    std::ofstream zyy((OutFilename + "_zyy.out").c_str());

    using hpx::lcos::unique_future;
    using hpx::async;
    using hpx::wait_all;
    std::vector<unique_future<jif3D::rvec> > ImplResult;
    ImplResult.reserve(nrealmax);
    Calc_action CalcImpl;
    std::vector<hpx::naming::id_type> localities = hpx::find_all_localities();

    for (size_t nreal = 0; nreal < nrealmax; ++nreal)
      {
        hpx::naming::id_type const locality_id = localities.at(nreal % localities.size());
        ImplResult.push_back(
            async(CalcImpl, locality_id, Model, bg_conductivity, phase1cond, phase2cond,
                phase1frac));
      }
    wait_all(ImplResult);

    for (size_t nreal = 0; nreal < nrealmax; ++nreal)
      {
        jif3D::rvec Impedances = ImplResult[nreal].get();
        std::string realstring(jif3D::stringify(nreal));
        zxx << nreal << " " << Impedances(0) << " " << Impedances(1) << std::endl;
        zxy << nreal << " " << Impedances(2) << " " << Impedances(3) << std::endl;
        zyx << nreal << " " << Impedances(4) << " " << Impedances(5) << std::endl;
        zyy << nreal << " " << Impedances(6) << " " << Impedances(7) << std::endl;

      }
    return hpx::finalize();
  }

int main(int argc, char* argv[])
  {
    return hpx::init(argc, argv);
  }
