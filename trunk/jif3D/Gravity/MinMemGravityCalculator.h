//============================================================================
// Name        : MinMemGravityCalculator.h
// Author      : Dec 10, 2008
// Version     :
// Copyright   : 2008, mmoorkamp
//============================================================================

#ifndef MINMEMGRAVITYCALCULATOR_H_
#define MINMEMGRAVITYCALCULATOR_H_

#include <boost/shared_ptr.hpp>
#include <boost/serialization/serialization.hpp>
#include <boost/serialization/base_object.hpp>
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
    private:
      friend class boost::serialization::access;
      //! Provide serialization to be able to store objects and, more importantly for simpler MPI parallelization
      template<class Archive>
      void serialize(Archive & ar, const unsigned int version)
        {
          ar & boost::serialization::base_object<ThreeDGravityCalculator>(*this);
        }
    public:
      //! The implementation of the forward calculation
      virtual rvec Calculate(const ThreeDGravityModel &Model);
      //! We have to implement this function even though it does not do anything
      virtual void HandleSensitivities(const size_t measindex)
        {
        }
      //! The constructor takes a shared pointer to an implementation object
      MinMemGravityCalculator(boost::shared_ptr<ThreeDGravityImplementation> TheImp);
      virtual ~MinMemGravityCalculator();
      };
  /* @} */
  }

#endif /* MINMEMGRAVITYCALCULATOR_H_ */
