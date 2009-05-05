//============================================================================
// Name        : FullSensitivityGravityCalculator.h
// Author      : Dec 12, 2008
// Version     :
// Copyright   : 2008, mmoorkamp
//============================================================================


#ifndef FULLSENSITIVITYGRAVITYCALCULATOR_H_
#define FULLSENSITIVITYGRAVITYCALCULATOR_H_

#include "CachedGravityCalculator.h"

namespace jiba
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
    class FullSensitivityGravityCalculator: public jiba::CachedGravityCalculator
      {
    private:
      rmat Sensitivities;
      //Calculates the raw gravity/FTG data from cached sensitivities without applying any transformation
      rvec CalculateRawData(const ThreeDGravityModel &Model);
      //Calculate a new model when no cached information is available or valid
      virtual rvec CalculateNewModel(const ThreeDGravityModel &Model);
      //calculate Data with applied transformation when cached sensitivities are still valid
      virtual rvec CalculateCachedResult(const ThreeDGravityModel &Model);
      virtual rvec CachedLQDerivative(const ThreeDGravityModel &Model, const rvec &Misfit);
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
      FullSensitivityGravityCalculator(boost::shared_ptr<
          ThreeDGravityImplementation> TheImp);
      virtual ~FullSensitivityGravityCalculator();
      };
  /* @} */
  }

#endif /* FULLSENSITIVITYGRAVITYCALCULATOR_H_ */
