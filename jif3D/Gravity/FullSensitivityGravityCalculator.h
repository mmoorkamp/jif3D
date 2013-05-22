//============================================================================
// Name        : FullSensitivityGravityCalculator.h
// Author      : Dec 12, 2008
// Version     :
// Copyright   : 2008, mmoorkamp
//============================================================================

#ifndef FULLSENSITIVITYGRAVITYCALCULATOR_H_
#define FULLSENSITIVITYGRAVITYCALCULATOR_H_

#include <boost/serialization/serialization.hpp>
#include <boost/serialization/base_object.hpp>
#include "CachedGravityCalculator.h"

namespace jif3D
  {
    /** \addtogroup gravity Gravity forward modeling, display and inversion */
    /* @{ */
    //! A calculator class that stores the full sensitivity matrix for inversion or fast consecutive forward calculations
    /*! This class stores the complete sensitivity matrix with all entries.
     * This allows to perform direct linear inversion and makes series of forward calculations
     * with identical geometries but changing density values extremely fast.
     * The disadvantage of this approach is that large models take up huge amounts of memory.
     * The actual management of the caching of calculations is performed by the base class
     * CachedGravityCalculator, this class only manages the storage and how cached results
     * are calculated.
     *
     * Important note: The stored sensitivities are the raw sensitivities irrespective of any
     * transformation. This is necessary to calculate correct gradients for transformed data.
     */
    class FullSensitivityGravityCalculator: public jif3D::CachedGravityCalculator
      {
    private:
      //! The full sensitivity matrix G for the gravity measurements, we use it to calculate the forward response through d = G*m
      rmat Sensitivities;
      //! Calculates the raw gravity/FTG data from cached sensitivities without applying any transformation
      virtual rvec CalculateRawData(const ThreeDGravityModel &Model);
      //! Calculates the raw gravity/FTG derivative from cached sensitivities without applying any transformation
      virtual rvec CalculateRawLQDerivative(const ThreeDGravityModel &Model,
          const rvec &Misfit);
      //! Calculate a new model when no cached information is available or valid
      virtual rvec CalculateNewModel(const ThreeDGravityModel &Model);
      //! calculate Data with applied transformation when cached sensitivities are still valid
      virtual rvec CalculateCachedResult(const ThreeDGravityModel &Model);
      //! calculate derivative of least squares objective function with applied transformation when cached sensitivities are still valid
      virtual rvec CachedLQDerivative(const ThreeDGravityModel &Model,
          const rvec &Misfit);
      friend class boost::serialization::access;
      //! Provide serialization to be able to store objects and, more importantly for simpler MPI parallelization
      template<class Archive>
      void serialize(Archive & ar, const unsigned int version)
        {
          ar & boost::serialization::base_object<ThreeDGravityCalculator>(*this);
          ar & Sensitivities;
        }
    public:
      //! return a read only copy of the sensitivity matrix, this guarantees that cache information is preserved
      const rmat &GetSensitivities() const
        {
          return Sensitivities;
        }
      //! For efficiency we sometimes operate directly on the sensitivities, as we don't have guaranteed cache information, we enforce recalculation
      rmat &SetSensitivities()
        {
          InvalidateCache();
          return Sensitivities;
        }
      //! This function is called by the implementation classes and allows to integrate the results from a single measurement
      virtual void HandleSensitivities(const size_t measindex);
      //! The constructor takes a shared pointer to an implementation object
      FullSensitivityGravityCalculator(
          boost::shared_ptr<ThreeDGravityImplementation> TheImp);
      virtual ~FullSensitivityGravityCalculator();
      };
  /* @} */
  }

#endif /* FULLSENSITIVITYGRAVITYCALCULATOR_H_ */
