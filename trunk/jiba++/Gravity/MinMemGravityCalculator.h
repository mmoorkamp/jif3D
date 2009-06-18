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
    /** \addtogroup gravity Gravity forward modeling, display and inversion */
    /* @{ */
    //! A calculator class that uses the minimum amount of memory, no caching, no sensitivity information
    /*! This is a calculator class that uses only the minimum amount of memory. It is suitable
     * when the forward response for a certain model geometry only has to be calculated once
     * or if the model is so big that the sensitivity matrix cannot be stored in memory or on disk any more.
     */
    class MinMemGravityCalculator: public jiba::ThreeDGravityCalculator
      {
    public:
      //! The implementation of the forward calculation
      virtual rvec Calculate(const ThreeDGravityModel &Model);
      //! We have to implement this function even though it does not do anything
      virtual void HandleSensitivities(const size_t measindex)
        {
        }
      ;
      //! The constructor takes a shared pointer to an implementation object
      MinMemGravityCalculator(
          boost::shared_ptr<ThreeDGravityImplementation> TheImp);
      virtual ~MinMemGravityCalculator();
      };
  /* @} */
  }

#endif /* MINMEMGRAVITYCALCULATOR_H_ */
