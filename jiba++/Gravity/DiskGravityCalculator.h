/*
 * DiskGravityCalculator.h
 *
 *  Created on: Jun 3, 2009
 *      Author: mmoorkamp
 */

#ifndef DISKGRAVITYCALCULATOR_H_
#define DISKGRAVITYCALCULATOR_H_

#include "FullSensitivityGravityCalculator.h"

namespace jiba
{
//! This class works similar to FullSensitivityCalculator, only that it stores the sensitivities on the disk, not in memory
class DiskGravityCalculator: public jiba::FullSensitivityGravityCalculator
{
private:
	std::string filename;
	virtual rvec CalculateRawData(const ThreeDGravityModel &Model);
	virtual rvec CalculateNewModel(const ThreeDGravityModel &Model);
	virtual rvec CalculateRawLQDerivative(const ThreeDGravityModel &Model, const rvec &Misfit);
public:
	virtual void HandleSensitivities(const size_t measindex);
	DiskGravityCalculator(boost::shared_ptr<ThreeDGravityImplementation> TheImp);
	virtual ~DiskGravityCalculator();
};

}

#endif /* DISKGRAVITYCALCULATOR_H_ */
