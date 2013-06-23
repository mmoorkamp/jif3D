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
#include "../Global/FileUtil.h"
#include "../Global/convert.h"
#include "../MT/X3DModel.h"
#include "../MT/X3DMTCalculator.h"
#include "../MT/ReadWriteImpedances.h"

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
	//set the cell sizes and allocate memory for the mesh
	Model.SetHorizontalCellSize(deltax, deltay, nx, ny);
	Model.SetZCellSizes().resize(boost::extents[nz]);
	fill_n(Model.SetZCellSizes().begin(), nz, deltaz);
	Model.SetConductivities().resize(boost::extents[nx][ny][nz]);
	//ask for a conductivity to fill the mesh with
	double bg_conductivity = 1.0;
	std::cout << "Background Conductivity: ";
	std::cin >> bg_conductivity;
	double phase1cond, phase2cond, phase1frac;
	std::cout << "Conductivity Phase 1:";
	std::cin >> phase1cond;
	std::cout << "Conductivity Phase 2:";
	std::cin >> phase2cond;
	std::cout << "Fraction Phase 1: ";
	std::cin >> phase1frac;
	srand48(time(0));
	double frequency;
	std::cout << "Frequency: ";
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
	std::vector<double> bg_conductivities(Model.GetZCellSizes().size(),
			bg_conductivity);
	Model.SetBackgroundThicknesses(bg_thicknesses);
	Model.SetBackgroundConductivities(bg_conductivities);
	size_t nrealmax;
	std::cout << "Realizations: ";
	std::cin >> nrealmax;

	std::ofstream zxy("zxy.out");
	std::ofstream zyx("zyx.out");
	std::ofstream zxx("zxx.out");
	std::ofstream zyy("zyy.out");
	for (size_t nreal = 0; nreal < nrealmax; ++nreal)
	{
		std::cout << "Realization: " << nreal << std::endl;
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
		std::string realstring(jif3D::stringify(nreal));
		Model.WriteNetCDF(OutFilename + realstring + ".nc");
		jif3D::X3DMTCalculator Calculator;
		jif3D::rvec Impedances(Calculator.Calculate(Model));
		jif3D::WriteImpedancesToNetCDF(OutFilename + realstring + "_data.nc",
				Model.GetFrequencies(), Model.GetMeasPosX(),
				Model.GetMeasPosY(), Model.GetMeasPosZ(), Impedances);
		jif3D::rvec Errors(Impedances.size(), 0.0);
		jif3D::WriteImpedancesToMtt(OutFilename + realstring + "_data.nc",
				Model.GetFrequencies(), Impedances, Errors);
		Model.WriteVTK(OutFilename + realstring + ".vtk");
		zxx << nreal << " " << Impedances(0) << " " << Impedances(1) << std::endl;
		zxy << nreal << " " << Impedances(2) << " " << Impedances(3) << std::endl;
		zyx << nreal << " " << Impedances(4) << " " << Impedances(5) << std::endl;
		zyy << nreal << " " << Impedances(6) << " " << Impedances(7) << std::endl;
	}

}
