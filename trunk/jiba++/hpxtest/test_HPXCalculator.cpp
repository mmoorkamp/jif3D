//============================================================================
// Name        : test_ReadWriteX3D.cpp
// Author      : Feb 17, 2009
// Version     :
// Copyright   : 2009, mmoorkamp
//============================================================================

#include <hpx/hpx_init.hpp>
#include "../Global/NumUtil.h"
#include "HPXMTCalculator.h"
#include "../MT/X3DModel.h"
#include "../MT/MTEquations.h"
#include "../MT/MT2DForward.h"
#include "../MT/ReadWriteImpedances.h"




int hpx_main(int argc, char* argv[])
{
	//create a 3D version of a 1D model
		const size_t xsize = 20;
		const size_t ysize = 20;
		const size_t zsize = 10;
		const size_t nbglayers = 5;
		jiba::X3DModel Model;

		Model.SetZCellSizes().resize(boost::extents[zsize]);

		Model.SetConductivities().resize(boost::extents[xsize][ysize][zsize]);
		std::vector<double> bg_thicknesses(nbglayers), bg_conductivities(nbglayers);

		const double deltax = 100.0;
		const double deltay = 100.0;
		const double deltaz = 100.0;
		const double freq = 1.0;
		const double cond = 0.01;
		Model.SetHorizontalCellSize(deltax, deltay, xsize, ysize);

		std::fill_n(Model.SetZCellSizes().origin(), zsize, deltaz);
		std::fill_n(Model.SetConductivities().origin(), xsize * ysize * zsize,
				cond);
		std::fill_n(bg_conductivities.begin(), nbglayers, cond);
		bg_conductivities.back() *= 1.001;
		std::fill_n(bg_thicknesses.begin(), nbglayers, 200.0);

		Model.SetBackgroundConductivities(bg_conductivities);
		Model.SetBackgroundThicknesses(bg_thicknesses);
		Model.SetFrequencies().push_back(freq);

		jiba::HPXMTCalculator Calculator;
		for (size_t i = 0; i < xsize/2; ++i)
			for (size_t j = 0; j < ysize/2; ++j)
			{
				Model.AddMeasurementPoint(Model.GetXCoordinates()[i] + deltax / 2.0,
						Model.GetYCoordinates()[j] + deltay / 2.0, 0.0);
			}

		jiba::rvec Impedance = Calculator.Calculate(Model);
		std::complex<double> HsImp = jiba::ImpedanceHalfspace(freq, cond);
		const double prec = 0.05;
		for (size_t i = 0; i < xsize * ysize; ++i)
		{
			std::cout << Impedance(i * 8 + 2) << " " << HsImp.real() << std::endl;
			std::cout << Impedance(i * 8 + 3) << " " << HsImp.imag() << std::endl;
			std::cout << Impedance(i * 8 + 4) << " " << -HsImp.real() << std::endl;
			std::cout << Impedance(i * 8 + 5) << " " << -HsImp.imag() << std::endl;

		}
    // Any HPX application logic goes here...
    return hpx::finalize();
}

int main(int argc, char* argv[])
{

    // Initialize HPX, run hpx_main as the first HPX thread, and
    // wait for hpx::finalize being called.
    return hpx::init(argc, argv);

}

