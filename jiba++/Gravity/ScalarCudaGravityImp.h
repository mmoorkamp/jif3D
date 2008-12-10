/*
 * ScalarCudaGravityImp.h
 *
 *  Created on: Dec 10, 2008
 *      Author: mmoorkamp
 */

#ifndef SCALARCUDAGRAVITYIMP_H_
#define SCALARCUDAGRAVITYIMP_H_

#include "ThreeDGravityImplementation.h"

namespace jiba {

class ScalarCudaGravityImp: public jiba::ThreeDGravityImplementation {
private:
	double *d_xcoord, *d_ycoord, *d_zcoord;
	double *d_xsize, *d_ysize, *d_zsize;
	double *d_result;
	static const size_t ndatapermeas = 1;
public:
	virtual size_t GetDataPerMeasurement() {
		return ndatapermeas;
	}
	virtual rvec CalcBackground(const double xmeas, const double ymeas,
			const double zmeas, const double xwidth, const double ywidth,
			const double zwidth, const ThreeDGravityModel &Model,
			rmat &Sensitivities);
	virtual rvec CalcGridded(const double x_meas, const double y_meas,
			const double z_meas, const ThreeDGravityModel &Model,
			rmat &Sensitivities);
	virtual rvec Calculate(const ThreeDGravityModel &Model,
			ThreeDGravityCalculator &Calculator);
	ScalarCudaGravityImp();
	virtual ~ScalarCudaGravityImp();
};

}

#endif /* SCALARCUDAGRAVITYIMP_H_ */
