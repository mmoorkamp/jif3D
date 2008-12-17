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
    class FullSensitivityGravityCalculator: public jiba::CachedGravityCalculator
      {
    private:
      rmat Sensitivities;
      virtual rvec CalculateNewModel(const ThreeDGravityModel &Model);
      virtual rvec CalculateCachedResult(const ThreeDGravityModel &Model);
    public:
      const rmat &GetSensitivities(){return Sensitivities;}
      virtual void HandleSensitivities(const size_t measindex);
      FullSensitivityGravityCalculator(boost::shared_ptr<ThreeDGravityImplementation> TheImp);
      virtual ~FullSensitivityGravityCalculator();
      };
    /* @} */
  }

#endif /* FULLSENSITIVITYGRAVITYCALCULATOR_H_ */
