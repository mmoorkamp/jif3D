//============================================================================
// Name        : MinMemGravityCalculator.h
// Author      : Dec 10, 2008
// Version     :
// Copyright   : 2008, mmoorkamp
//============================================================================


#ifndef MINMEMGRAVITYCALCULATOR_H_
#define MINMEMGRAVITYCALCULATOR_H_

#include "ThreeDGravityCalculator.h"

namespace jiba
  {

    class MinMemGravityCalculator: public jiba::ThreeDGravityCalculator
      {
    public:

      virtual rvec Calculate(const ThreeDGravityModel &Model);
      MinMemGravityCalculator(ThreeDGravityImplementation &TheImp);

      virtual void HandleSensitivities(const size_t measindex){};
      MinMemGravityCalculator();
      virtual ~MinMemGravityCalculator();
      };

  }

#endif /* MINMEMGRAVITYCALCULATOR_H_ */
