//============================================================================
// Name        : EqualGeometry.h
// Author      : Mar 3, 2010
// Version     : 
// Copyright   : 2010, mmoorkamp
//============================================================================

#ifndef EQUALGEOMETRY_H_
#define EQUALGEOMETRY_H_

#include <iostream>
#include "../Global/NumUtil.h"
#include "../Global/FatalException.h"
#include "ThreeDModelBase.h"

namespace jiba {
//! Test if two models have the same grid geometry
/*! This function can be used to test whether two 3D models have the same grid geometry, i.e. the same number of cells
 * in each direction with the same sizes.
 * @param Model1 The first model for the comparison
 * @param Model2 the second model for the comparison
 * @return True if the geometries are equal, false otherwise
 */
inline bool EqualGridGeometry(const jiba::ThreeDModelBase &Model1,
		const jiba::ThreeDModelBase &Model2) {
	//we first make a cheap test and then the more expensive ones
	if (Model1.GetNModelElements() != Model2.GetNModelElements())
		return false;
	// as soon as one fails we return false
	//we check the sizes of the cells in all three directions
	const std::size_t nx = Model1.GetXCellSizes().num_elements();

	for (std::size_t i = 0; i < nx; ++i) {
		if (!roughlyEqual(Model1.GetXCellSizes()[i],
				Model2.GetXCellSizes()[i])) {
			std::cerr << "Cell sizes in x-direction do not match "
					<< Model1.GetXCellSizes()[i] << " "
					<< Model1.GetXCellSizes()[i] << std::endl;
			return false;
		}
	}
	for (std::size_t i = 0; i < nx; ++i) {
		if (!roughlyEqual(Model1.GetYCellSizes()[i],
				Model2.GetYCellSizes()[i])) {
			std::cerr << "Cell sizes in y-direction do not match "
					<< Model1.GetYCellSizes()[i] << " "
					<< Model1.GetYCellSizes()[i] << std::endl;
			return false;
		}
	}
	for (std::size_t i = 0; i < nx; ++i) {
		if (!roughlyEqual(Model1.GetZCellSizes()[i],
				Model2.GetZCellSizes()[i])) {
			std::cerr << "Cell sizes in z-direction do not match "
					<< Model1.GetZCellSizes()[i] << " "
					<< Model1.GetZCellSizes()[i] << std::endl;
			return false;
		}
	}
	//we only get here if everything is equal
	return true;
}
}
#endif /* EQUALGEOMETRY_H_ */
