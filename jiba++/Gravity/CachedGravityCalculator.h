//============================================================================
// Name        : CachedGravityCalculator.h
// Author      : Dec 12, 2008
// Version     :
// Copyright   : 2008, mmoorkamp
//============================================================================


#ifndef CACHEDGRAVITYCALCULATOR_H_
#define CACHEDGRAVITYCALCULATOR_H_

#include "ThreeDGravityCalculator.h"

namespace jiba
  {
    /** \addtogroup gravity Gravity forward modelling, display and inversion */
      /* @{ */
    class CachedGravityCalculator: public jiba::ThreeDGravityCalculator
      {
    private:
      bool HaveCache;
      ThreeDGravityModel::t3DModelDim OldXSizes;
      ThreeDGravityModel::t3DModelDim OldYSizes;
      ThreeDGravityModel::t3DModelDim OldZSizes;
      ThreeDGravityModel::tMeasPosVec OldMeasPosX;
      ThreeDGravityModel::tMeasPosVec OldMeasPosY;
      ThreeDGravityModel::tMeasPosVec OldMeasPosZ;
      void CopySizes(const ThreeDGravityModel::t3DModelDim &NewXSizes,
          const ThreeDGravityModel::t3DModelDim &NewYSizes,
          const ThreeDGravityModel::t3DModelDim &NewZSizes);
      void CopyMeasPos(const ThreeDGravityModel::tMeasPosVec &NewMeasPosX,
          const ThreeDGravityModel::tMeasPosVec &NewMeasPosY,
          const ThreeDGravityModel::tMeasPosVec &NewMeasPosZ);
      bool CheckGeometryChange(const ThreeDGravityModel &Model);
      bool CheckMeasPosChange(const ThreeDGravityModel &Model);
      virtual rvec CalculateNewModel(const ThreeDGravityModel &Model) = 0;
      virtual rvec CalculateCachedResult(const ThreeDGravityModel &Model) = 0;
    public:
      virtual rvec Calculate(const ThreeDGravityModel &Model);
      CachedGravityCalculator(
          boost::shared_ptr<ThreeDGravityImplementation> TheImp);
      virtual ~CachedGravityCalculator();
      };
    /* @} */
  }

#endif /* CACHEDGRAVITYCALCULATOR_H_ */
