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
#include "ThreeDGravMagCalculator.h"

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
    template<class ThreeDModelType>
    class CachedGravMagCalculator: public jif3D::ThreeDGravMagCalculator<ThreeDModelType>
      {
    private:
      //! Have we performed a calculation before and build up cached information
      bool HaveCache;
      //we have to store the model and measurement information of the previous
      //calculation, so we can decide whether the geometries have changed
      typedef typename ThreeDModelType::t3DModelDim ModelDimType;
      typedef typename ThreeDModelType::tMeasPosVec MeasPosType;
      typedef typename ThreeDModelType::tBackgroundVec BackgroundType;
      typename ThreeDModelType::ModelCacheType ModelCache;
      friend class boost::serialization::access;
      //! Provide serialization to be able to store objects and, more importantly for simpler MPI parallelization
      template<class Archive>
      void serialize(Archive & ar, const unsigned int version)
        {
          ar & boost::serialization::base_object<ThreeDGravMagCalculator>(*this);
        }
      //! The function declaration for the calculation of a new sensitivity matrix and rebuilding of caching information
      virtual rvec CalculateNewModel(const ThreeDModelType &Model) = 0;
      //! The function declaration for the calculation of gravity data using the cached information
      virtual rvec CalculateCachedResult(const ThreeDModelType &Model) = 0;
      //! The function declaration for the calculation of the gradient using the cached information
      virtual rvec CachedLQDerivative(const ThreeDModelType &Model,
          const rvec &Misfit) = 0;
    protected:
      //! Through this function derived classes can signal that the sensitivity information is not valid any more
      void InvalidateCache()
        {
          HaveCache = false;
        }
    public:
      //! Calculate the data for the given gravity model
      virtual rvec Calculate(const ThreeDModelType &Model);
      //! We overwrite the base class method to use the caching information
      virtual rvec LQDerivative(const ThreeDModelType &Model, const rvec &Misfit);
      //! The constructor takes a shared pointer to an implementation object
      CachedGravMagCalculator(
          boost::shared_ptr<ThreeDGravMagImplementation<ThreeDModelType> > TheImp);
      virtual ~CachedGravMagCalculator();
      };

    template<class ThreeDModelType>
    CachedGravMagCalculator<ThreeDModelType>::CachedGravMagCalculator(
        boost::shared_ptr<ThreeDGravMagImplementation<ThreeDModelType> > TheImp) :
        ThreeDGravMagCalculator<ThreeDModelType>(TheImp), HaveCache(false)
      {

      }

    template<class ThreeDModelType>
    CachedGravMagCalculator<ThreeDModelType>::~CachedGravMagCalculator()
      {
      }



    /*! This function compares the geometry and measurements of the passed
     * model with the information stored during the last call and either
     * calls CalculateCachedResult or CalculateNewModel which have to
     * be implemented in a derived class.
     * @param Model The 3D gravity model
     * @return The forward response for this model
     */
    template<class ThreeDModelType>
    rvec CachedGravMagCalculator<ThreeDModelType>::Calculate(const ThreeDModelType &Model)
      {
        //check that all modeling information is consistent
        CheckModelConsistency(Model);
        //check whether the model geometry has changed since the last calculation
        bool HasChanged = ModelCache.HasChanged(Model);
        //if we have cached information and nothing has changed, we use the cache
        if (HaveCache && !HasChanged)
          {
            return CalculateCachedResult(Model);
          }
        //we only get here if we need to recalculate
        //we have to make sure the calculation finishes properly before we can guarantee the cache and return the result
        ThreeDGravMagCalculator<ThreeDModelType>::SetCurrentSensitivities().resize(
            ThreeDGravMagCalculator<ThreeDModelType>::Imp.get()->RawDataPerMeasurement(),
            Model.GetDensities().num_elements() + Model.GetBackgroundDensities().size());
        rvec result = CalculateNewModel(Model);
        HaveCache = true;

        return result;
      }

    template<class ThreeDModelType>
    rvec CachedGravMagCalculator<ThreeDModelType>::LQDerivative(
        const ThreeDModelType &Model, const rvec &Misfit)
      {
        //check that all modeling information is consistent
        CheckModelConsistency(Model);
        //check whether the model geometry has changed since the last calculation
        bool HasChanged = ModelCache.HasChanged(Model);
        //if we have cached information and nothing has changed, we use the cache
        if (HaveCache && !HasChanged)
          {
            return CachedLQDerivative(Model, Misfit);
          }
        //we only get here if we need to recalculate
        //we have to make sure the calculation finishes properly before we can guarantee the cache and return the result
        ThreeDGravMagCalculator<ThreeDModelType>::SetCurrentSensitivities().resize(
            ThreeDGravMagCalculator<ThreeDModelType>::Imp.get()->GetDataPerMeasurement(),
            Model.GetDensities().num_elements() + Model.GetBackgroundDensities().size());
        return ThreeDGravMagCalculator<ThreeDModelType>::LQDerivative(Model, Misfit);

      }
  /* @} */
  }

#endif /* CACHEDGRAVITYCALCULATOR_H_ */
