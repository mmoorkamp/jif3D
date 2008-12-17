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
    /** \addtogroup gravity Gravity forward modelling, display and inversion */
      /* @{ */
    //! A calculator class that uses the minimum amount of memory, no caching, no sensitivity information
    class MinMemGravityCalculator: public jiba::ThreeDGravityCalculator
      {
    public:
      //! The implementation of the forward calculation
      virtual rvec Calculate(const ThreeDGravityModel &Model);
      MinMemGravityCalculator(boost::shared_ptr<ThreeDGravityImplementation> TheImp);
      //! We have to implement this function even though it does not do anything
      virtual void HandleSensitivities(const size_t measindex){};
      MinMemGravityCalculator();
      virtual ~MinMemGravityCalculator();
      };
    /* @} */
  }

#endif /* MINMEMGRAVITYCALCULATOR_H_ */
