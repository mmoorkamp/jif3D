/*
 * ScalarCudaGravityImp.cpp
 *
 *  Created on: Dec 10, 2008
 *      Author: mmoorkamp
 */

#include "ScalarCudaGravityImp.h"
#include <numeric>
#include "ThreeDGravityCalculator.h"
namespace jiba {

extern "C"void SingleScalarMeas(const double x_meas, const double y_meas,
		const double z_meas, double *d_xcoord, double *d_ycoord,
		double *d_zcoord, double *d_xsize, double *d_ysize, double *d_zsize,
		double *d_result, const int nx, const int ny, const int nz,
		double *returnvalue);

extern "C" void PrepareData(double *d_xcoord, double *d_ycoord, double *d_zcoord,
		double *d_xsize, double *d_ysize, double *d_zsize, double *d_result,
		const double *xcoord,const double *ycoord,const double *zcoord,
		const double *xsize,const double *ysize,const double *zsize, unsigned int nx,unsigned int ny,unsigned int nz);

extern "C" void FreeData(double *d_xcoord, double *d_ycoord, double *d_zcoord,
		double *d_xsize, double *d_ysize, double *d_zsize, double *d_result);

ScalarCudaGravityImp::ScalarCudaGravityImp() {
	// TODO Auto-generated constructor stub

}

ScalarCudaGravityImp::~ScalarCudaGravityImp() {
	// TODO Auto-generated destructor stub
}

rvec ScalarCudaGravityImp::CalcBackground(const double xmeas, const double ymeas,
		const double zmeas, const double xwidth, const double ywidth,
		const double zwidth, const ThreeDGravityModel &Model,
		rmat &Sensitivities)
{

}

rvec ScalarCudaGravityImp::CalcGridded(const double x_meas, const double y_meas,
		const double z_meas, const ThreeDGravityModel &Model,
		rmat &Sensitivities)
{
	SingleScalarMeas(x_meas,y_meas,z_meas,d_xcoord,d_ycoord,d_zcoord,d_xsize,d_ysize,d_zsize,d_result,
			Model.GetDensities().shape()[0],Model.GetDensities().shape()[1],Model.GetDensities().shape()[2],&Sensitivities.data()[0]);
	rvec result(1);
	result(0) = std::inner_product(Sensitivities.data().begin(),
			Sensitivities.data().end(),Model.GetDensities().origin(),0.0);
	return result;
}

rvec ScalarCudaGravityImp::Calculate(const ThreeDGravityModel &Model,ThreeDGravityCalculator &Calculator)
{

	const unsigned int nx = Model.GetDensities().shape()[0];
	const unsigned int ny = Model.GetDensities().shape()[1];
	const unsigned int nz = Model.GetDensities().shape()[2];
	const unsigned int nelements = nx * ny * nz;
	PrepareData(d_xcoord,d_ycoord,d_zcoord,d_xsize,d_ysize,d_zsize,d_result,Model.GetXCoordinates().data(),Model.GetYCoordinates().data(),
			Model.GetZCoordinates().data(),Model.GetXCellSizes().data(),Model.GetYCellSizes().data(),Model.GetZCellSizes().data(),nx,ny,nz);
	const unsigned int nmeas = Model.GetMeasPosX().size();
	rvec result(nmeas);
	for (size_t i = 0; i < nmeas; ++i)
	{
		rvec currresult = CalcGridded(Model.GetMeasPosX()[i],Model.GetMeasPosY()[i],Model.GetMeasPosZ()[i],Model,Calculator.SetCurrentSensitivities());

		result(i) = currresult(0);

	}
	FreeData(d_xcoord,d_ycoord,d_zcoord,d_xsize,d_ysize,d_zsize,d_result);
	return result;
}

}
