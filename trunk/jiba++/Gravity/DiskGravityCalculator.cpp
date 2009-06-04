/*
 * DiskGravityCalculator.cpp
 *
 *  Created on: Jun 3, 2009
 *      Author: mmoorkamp
 */

#include <unistd.h>
#include <boost/filesystem.hpp>
#include "../Global/convert.h"
#include "../Global/VecMat.h"
#include "DiskGravityCalculator.h"
#include <fstream>

namespace jiba
{

rvec DiskGravityCalculator::CalculateNewModel(
    const ThreeDGravityModel &Model)
  {
	boost::filesystem::remove_all(filename);
    //then forward the call to the implementation object
    return Imp.get()->Calculate(Model, *this);
  }

void DiskGravityCalculator::HandleSensitivities(const size_t measindex)
{
	std::fstream outfile(filename.c_str(), std::ios::out | std::ios::binary
			| std::ios::app);
	jiba::rmat CurrSens(SetCurrentSensitivities());

	for (size_t i = 0; i < CurrSens.size1(); ++i)
	{
	jiba::rvec Row(ublas::matrix_row<jiba::rmat>(CurrSens,i));
	outfile.write(reinterpret_cast<char *> (&Row[0]), Row.size() * sizeof(double));
	}
}

rvec DiskGravityCalculator::CalculateRawData(const ThreeDGravityModel &Model)
{
	std::fstream infile(filename.c_str(),std::ios::in|std::ios::binary);

    const size_t nmeas = Model.GetMeasPosX().size()
        * Imp.get()->RawDataPerMeasurement();
    const size_t ngrid = Model.GetDensities().num_elements();
    const size_t nmod = ngrid + Model.GetBackgroundThicknesses().size();

    rvec DensVector(nmod);
    //copy the 3D model structure and the background into a vector of densities
    std::copy(Model.GetDensities().origin(), Model.GetDensities().origin()
        + ngrid, DensVector.begin());
    std::copy(Model.GetBackgroundDensities().begin(),
        Model.GetBackgroundDensities().end(), DensVector.begin() + ngrid);
    rvec result(nmeas);
    rvec CurrSens(nmod);
    for (size_t i = 0; i < nmeas; ++i)
    {
    	infile.read(reinterpret_cast<char *> (&CurrSens.data()[0]),nmod*sizeof(double));
    	result(i) = boost::numeric::ublas::inner_prod(CurrSens,DensVector);
    }
    return result;

}

rvec DiskGravityCalculator::CalculateRawLQDerivative(const ThreeDGravityModel &Model, const rvec &Misfit)
{
	std::fstream infile(filename.c_str(),std::ios::in|std::ios::binary);

	    const size_t nmeas = Model.GetMeasPosX().size()
	        * Imp.get()->RawDataPerMeasurement();
	    const size_t ngrid = Model.GetDensities().num_elements();
	    const size_t nmod = ngrid + Model.GetBackgroundThicknesses().size();
	    rvec result(nmod);
	    result.clear();
	    rvec CurrSens(nmod);
	    for (size_t i = 0; i < nmeas; ++i)
	    {
	    	infile.read(reinterpret_cast<char *> (&CurrSens.data()[0]),nmod*sizeof(double));
	    	result  += CurrSens * Misfit(i);
	    }
	    return result;
}

DiskGravityCalculator::DiskGravityCalculator(boost::shared_ptr<ThreeDGravityImplementation> TheImp) :
   FullSensitivityGravityCalculator(TheImp)
{
	filename = "grav" + jiba::stringify(getpid()) + jiba::stringify(this);
}

DiskGravityCalculator::~DiskGravityCalculator()
{
	boost::filesystem::remove_all(filename);

}

}
