//============================================================================
// Name        : OMPMagneticImp.cpp
// Author      : 6 Nov 2013
// Version     : 
// Copyright   : 2013, mm489
//============================================================================

#include "OMPMagneticImp.h"
#include "../Gravity/BasicGravElements.h"

namespace jif3D
{

rvec OMPMagneticImp::CalcGridded(const size_t measindex,
		const ThreeDMagneticModel &Model, rmat &Sensitivities)
{
	const double BxComp = cos(Inclination) * cos(Declination);
	const double ByComp = cos(Inclination) * sin(Declination);
	const double BzComp = sin(Inclination);

	const size_t xsize = Model.GetSusceptibilities().shape()[0];
	const size_t ysize = Model.GetSusceptibilities().shape()[1];
	const size_t zsize = Model.GetSusceptibilities().shape()[2];
	const double x_meas = Model.GetMeasPosX()[measindex];
	const double y_meas = Model.GetMeasPosY()[measindex];
	const double z_meas = Model.GetMeasPosZ()[measindex];
	const int nmod = xsize * ysize * zsize;
	const bool storesens = (Sensitivities.size1() >= ndatapermeas)
			&& (Sensitivities.size2() >= size_t(nmod));

	rmat currvalue(3, 3);

	//we cannot add up a user defined quantity in parallel
	//so break up the tensor into its component with different variables
	//and assign the results after the parallel loop
	double Bx = 0.0, By = 0.0, Bz = 0.0;

	//sum up the contributions of all prisms
#pragma omp parallel default(shared) private(currvalue) reduction(+:Bx,By,Bz)
	{
		//instead of nested loops over each dimension, we have one big
		//loop over all elements, this allows for a nearly infinite number
		//of parallel processors
#pragma omp for
		for (int offset = 0; offset < nmod; ++offset)
		{
			int xindex, yindex, zindex;
			//we still need the indices for each dimension
			//so we have to convert our loop variable
			Model.OffsetToIndex(offset, xindex, yindex, zindex);
			// we reuse the calculation for the FTG matrix, as the equations are
			//identical
			//currvalue contains only the geometric term
			currvalue = CalcTensorBoxTerm(x_meas, y_meas, z_meas,
					XCoord[xindex], YCoord[yindex], ZCoord[zindex],
					XSizes[xindex], YSizes[yindex], ZSizes[zindex]);
			//to we have to multiply each element by the density
			const double Susceptibility =
					Model.GetSusceptibilities()[xindex][yindex][zindex];
			const double BxSens = (currvalue(0, 0) * BxComp
					+ currvalue(0, 1) * ByComp + currvalue(0, 2) * BzComp)
					* FieldStrength;
			const double BySens = (currvalue(1, 0) * BxComp
					+ currvalue(1, 1) * ByComp + currvalue(2, 2) * BzComp)
					* FieldStrength;
			const double BzSens = (currvalue(2, 0) * BxComp
					+ currvalue(2, 1) * ByComp + currvalue(2, 2) * BzComp)
					* FieldStrength;
			Bx += BxSens * Susceptibility;
			By += BySens * Susceptibility;
			Bz += BzSens * Susceptibility;
			if (storesens)
			{
				Sensitivities(0, offset) = BxSens;
				Sensitivities(1, offset) = BySens;
				Sensitivities(2, offset) = BzSens;
			}
		}
	} //end of parallel region

	rvec returnvalue(ndatapermeas);
	returnvalue(0) = Bx;
	returnvalue(1) = By;
	returnvalue(2) = Bz;

	return returnvalue;

}

OMPMagneticImp::OMPMagneticImp(double Inc, double Dec, double Fs) :
		Inclination(Inc), Declination(Dec), FieldStrength(Fs)
{
	// TODO Auto-generated constructor stub

}

OMPMagneticImp::~OMPMagneticImp()
{
	// TODO Auto-generated destructor stub
}

} /* namespace jif3D */
