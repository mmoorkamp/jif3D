//============================================================================
// Name        : ThreeDGravityImplementation.h
// Author      : Dec 10, 2008
// Version     :
// Copyright   : 2008, mmoorkamp
//============================================================================


#ifndef THREEDGRAVITYIMPLEMENTATION_H_
#define THREEDGRAVITYIMPLEMENTATION_H_

#include "ThreeDGravityModel.h"
#include "../Global/VecMat.h"

namespace jiba
  {
	class ThreeDGravityCalculator;
    class ThreeDGravityImplementation
      {
    public:
      virtual size_t GetDataPerMeasurement() = 0;
      //! Calculate the response of the 1D background
      virtual rvec CalcBackground(const double xmeas, const double ymeas,
          const double zmeas, const double xwidth, const double ywidth,
          const double zwidth, const ThreeDGravityModel &Model,
          rmat &Sensitivities) = 0;
      virtual rvec CalcGridded(const double x_meas, const double y_meas,
          const double z_meas, const ThreeDGravityModel &Model,
          rmat &Sensitivities) = 0;
      virtual rvec Calculate(const ThreeDGravityModel &Model,ThreeDGravityCalculator &Calculator);
      ThreeDGravityImplementation();
      virtual ~ThreeDGravityImplementation();
      };

  }
#endif /* THREEDGRAVITYIMPLEMENTATION_H_ */
