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
    /** \addtogroup gravity Gravity forward modeling, display and inversion */
    /* @{ */
    //! The base class for all calculator classes that use sensitivity information to accelerate consecutive model calculations
    /*! This class analyzes the geometry of the model and measurements each time Calculate() is called.
     * If the geometries have not changed since the last call, CalculateCachedResult() is called where a derived
     * class can implement an accelerated forward calculation using information acquired during the previous
     * calculation. Otherwise CalculateNewModel() is called and the derived class has to calculate a new
     * model and rebuild the caching information.
     *
     * Examples for derived classes are FullSensitivityGravityCalculator and WaveletCompressedGravityCalculator.
     */
    class CachedGravityCalculator: public jiba::ThreeDGravityCalculator
      {
    private:
      //have we performed a calculation before and build up cached information
      bool HaveCache;
      //we have to store the model and measurement information of the previous
      //calculation, so we can decide whether the geometries have changed
      ThreeDGravityModel::t3DModelDim OldXSizes;
      ThreeDGravityModel::t3DModelDim OldYSizes;
      ThreeDGravityModel::t3DModelDim OldZSizes;
      ThreeDGravityModel::tMeasPosVec OldMeasPosX;
      ThreeDGravityModel::tMeasPosVec OldMeasPosY;
      ThreeDGravityModel::tMeasPosVec OldMeasPosZ;
      ThreeDGravityModel::tBackgroundVec OldBackgroundThick;
      ThreeDGravityModel::tBackgroundVec OldBackgroundDens;
      void CopySizes(const ThreeDGravityModel::t3DModelDim &NewXSizes,
          const ThreeDGravityModel::t3DModelDim &NewYSizes,
          const ThreeDGravityModel::t3DModelDim &NewZSizes);
      void CopyMeasPos(const ThreeDGravityModel::tMeasPosVec &NewMeasPosX,
          const ThreeDGravityModel::tMeasPosVec &NewMeasPosY,
          const ThreeDGravityModel::tMeasPosVec &NewMeasPosZ);
      bool CheckGeometryChange(const ThreeDGravityModel &Model);
      bool CheckMeasPosChange(const ThreeDGravityModel &Model);
      bool CheckBackgroundChange(const ThreeDGravityModel &Model);
      //! The function declaration for the calculation of a new sensitivity matrix and rebuilding of caching information
      virtual rvec CalculateNewModel(const ThreeDGravityModel &Model) = 0;
      //! The function declaration for the calculation of gravity data using the cached information
      virtual rvec CalculateCachedResult(const ThreeDGravityModel &Model) = 0;
    protected:
      //! Through this function derived classes can signal that the sensitivity information is not valid any more
      void InvalidateCache()
        {
          HaveCache = false;
        }
    public:
      //! Calculate the data for the given gravity model
      virtual rvec Calculate(const ThreeDGravityModel &Model);
      //! The constructor takes a shared pointer to an implementation object
      CachedGravityCalculator(
          boost::shared_ptr<ThreeDGravityImplementation> TheImp);
      virtual ~CachedGravityCalculator();
      };
  /* @} */
  }

#endif /* CACHEDGRAVITYCALCULATOR_H_ */
