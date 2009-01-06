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
    /** \addtogroup gravity Gravity forward modelling, display and inversion */
      /* @{ */
    //! A calculator class that stores the full sensitivity matrix for inversion or fast consecutive forward calculations
    /*! This class stores the complete sensitivity matrix with all entries.
     * This allows to perform direct linear inversion and makes series of forward calculations
     * with identical geometries but changing density values extremely fast.
     * The disadvantage of this approach is that large models take up huge amounts of memory.
     * The actual management of the caching of calculations is performed by the base class
     * CachedGravityCalculator, this class only manages the storage and how cached results
     * are calculated.
     */
    class FullSensitivityGravityCalculator: public jiba::CachedGravityCalculator
      {
    private:
      rmat Sensitivities;
      virtual rvec CalculateNewModel(const ThreeDGravityModel &Model);
      virtual rvec CalculateCachedResult(const ThreeDGravityModel &Model);
    public:
      //! return a read only copy of the sensitivity matrix
      const rmat &GetSensitivities(){return Sensitivities;}
      //! This function is called by the implementation classes and allows to integrate the results from a single measurement
      virtual void HandleSensitivities(const size_t measindex);
      FullSensitivityGravityCalculator(boost::shared_ptr<ThreeDGravityImplementation> TheImp);
      virtual ~FullSensitivityGravityCalculator();
      };
    /* @} */
  }

#endif /* FULLSENSITIVITYGRAVITYCALCULATOR_H_ */
