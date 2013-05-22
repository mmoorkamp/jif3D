//============================================================================
// Name        : CachedGravityCalculator.h
// Author      : Dec 12, 2008
// Version     :
// Copyright   : 2008, mmoorkamp
//============================================================================


#ifndef CACHEDGRAVITYCALCULATOR_H_
#define CACHEDGRAVITYCALCULATOR_H_

#include <boost/serialization/serialization.hpp>
#include <boost/serialization/base_object.hpp>
#include "ThreeDGravityCalculator.h"

namespace jif3D
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
    class CachedGravityCalculator: public jif3D::ThreeDGravityCalculator
      {
    private:
      //! Have we performed a calculation before and build up cached information
      bool HaveCache;
      //we have to store the model and measurement information of the previous
      //calculation, so we can decide whether the geometries have changed
      //! The size of the model cells in x-direction for the previous calculations
      ThreeDGravityModel::t3DModelDim OldXSizes;
      //! The size of the model cells in y-direction for the previous calculations
      ThreeDGravityModel::t3DModelDim OldYSizes;
      //! The size of the model cells in z-direction for the previous calculations
      ThreeDGravityModel::t3DModelDim OldZSizes;
      //! The measurement coordinates in x-direction for the previous calculations
      ThreeDGravityModel::tMeasPosVec OldMeasPosX;
      //! The measurement coordinates in y-direction for the previous calculations
      ThreeDGravityModel::tMeasPosVec OldMeasPosY;
      //! The measurement coordinates in z-direction for the previous calculations
      ThreeDGravityModel::tMeasPosVec OldMeasPosZ;
      //! The thickness of the background layers for the previous calculation
      ThreeDGravityModel::tBackgroundVec OldBackgroundThick;
      //! The density of the background layers for the previous calculation
      ThreeDGravityModel::tBackgroundVec OldBackgroundDens;
      friend class boost::serialization::access;
      //! Provide serialization to be able to store objects and, more importantly for simpler MPI parallelization
      template<class Archive>
      void serialize(Archive & ar, const unsigned int version)
        {
          ar & boost::serialization::base_object<ThreeDGravityCalculator>(*this);
          ar & OldXSizes;
          ar & OldYSizes;
          ar & OldZSizes;
          ar & OldMeasPosX;
          ar & OldMeasPosY;
          ar & OldMeasPosZ;
          ar & OldBackgroundThick;
          ar & OldBackgroundDens;
        }
      //! Copy the cell sizes from the current model to store them for caching
      void CopySizes(const ThreeDGravityModel::t3DModelDim &NewXSizes,
          const ThreeDGravityModel::t3DModelDim &NewYSizes,
          const ThreeDGravityModel::t3DModelDim &NewZSizes);
      //! Copy the measurement positions from the current model to store them for caching
      void CopyMeasPos(const ThreeDGravityModel::tMeasPosVec &NewMeasPosX,
          const ThreeDGravityModel::tMeasPosVec &NewMeasPosY,
          const ThreeDGravityModel::tMeasPosVec &NewMeasPosZ);
      //! Check whether the model geometry, i.e. cell sizes, has changed since the last calculation
      bool CheckGeometryChange(const ThreeDGravityModel &Model);
      //! Check whether the measurement positions have changed since the last calculation
      bool CheckMeasPosChange(const ThreeDGravityModel &Model);
      //! Check wether the 1D background has changed since the last calculation
      bool CheckBackgroundChange(const ThreeDGravityModel &Model);
      //! The function declaration for the calculation of a new sensitivity matrix and rebuilding of caching information
      virtual rvec CalculateNewModel(const ThreeDGravityModel &Model) = 0;
      //! The function declaration for the calculation of gravity data using the cached information
      virtual rvec CalculateCachedResult(const ThreeDGravityModel &Model) = 0;
      //! The function declaration for the calculation of the gradient using the cached information
      virtual rvec CachedLQDerivative(const ThreeDGravityModel &Model, const rvec &Misfit) = 0;
    protected:
      //! Through this function derived classes can signal that the sensitivity information is not valid any more
      void InvalidateCache()
        {
          HaveCache = false;
        }
    public:
      //! Calculate the data for the given gravity model
      virtual rvec Calculate(const ThreeDGravityModel &Model);
      //! We overwrite the base class method to use the caching information
      virtual rvec LQDerivative(const ThreeDGravityModel &Model, const rvec &Misfit);
      //! The constructor takes a shared pointer to an implementation object
      CachedGravityCalculator(
          boost::shared_ptr<ThreeDGravityImplementation> TheImp);
      virtual ~CachedGravityCalculator();
      };
  /* @} */
  }

#endif /* CACHEDGRAVITYCALCULATOR_H_ */
