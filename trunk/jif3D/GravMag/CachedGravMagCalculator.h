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
      //! The size of the model cells in x-direction for the previous calculations
      ModelDimType OldXSizes;
      //! The size of the model cells in y-direction for the previous calculations
      ModelDimType OldYSizes;
      //! The size of the model cells in z-direction for the previous calculations
      ModelDimType OldZSizes;
      //! The measurement coordinates in x-direction for the previous calculations
      MeasPosType OldMeasPosX;
      //! The measurement coordinates in y-direction for the previous calculations
      MeasPosType OldMeasPosY;
      //! The measurement coordinates in z-direction for the previous calculations
      MeasPosType OldMeasPosZ;
      //! The thickness of the background layers for the previous calculation
      BackgroundType OldBackgroundThick;
      //! The density of the background layers for the previous calculation
      BackgroundType OldBackgroundDens;
      friend class boost::serialization::access;
      //! Provide serialization to be able to store objects and, more importantly for simpler MPI parallelization
      template<class Archive>
      void serialize(Archive & ar, const unsigned int version)
        {
          ar & boost::serialization::base_object<ThreeDGravMagCalculator>(*this);
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
      void CopySizes(const ModelDimType &NewXSizes, const ModelDimType &NewYSizes,
          const ModelDimType &NewZSizes);
      //! Copy the measurement positions from the current model to store them for caching
      void CopyMeasPos(const MeasPosType &NewMeasPosX, const MeasPosType &NewMeasPosY,
          const MeasPosType &NewMeasPosZ);
      //! Check whether the model geometry, i.e. cell sizes, has changed since the last calculation
      bool CheckGeometryChange(const ThreeDModelType &Model);
      //! Check whether the measurement positions have changed since the last calculation
      bool CheckMeasPosChange(const ThreeDModelType &Model);
      //! Check wether the 1D background has changed since the last calculation
      bool CheckBackgroundChange(const ThreeDModelType &Model);
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
        ThreeDGravMagCalculator<ThreeDModelType>(TheImp), HaveCache(false), OldXSizes(), OldYSizes(), OldZSizes(), OldMeasPosX(), OldMeasPosY(), OldMeasPosZ(), OldBackgroundThick(), OldBackgroundDens()
      {

      }

    template<class ThreeDModelType>
    CachedGravMagCalculator<ThreeDModelType>::~CachedGravMagCalculator()
      {
      }

    template<class ThreeDModelType>
    void CachedGravMagCalculator<ThreeDModelType>::CopySizes(
        const ModelDimType &NewXSizes, const ModelDimType &NewYSizes,
        const ModelDimType &NewZSizes)
      {
        //make sure we have enough memory
        OldXSizes.resize(boost::extents[NewXSizes.size()]);
        OldYSizes.resize(boost::extents[NewYSizes.size()]);
        OldZSizes.resize(boost::extents[NewZSizes.size()]);
        //copy old sizes into member variables for next comparison
        std::copy(NewXSizes.begin(), NewXSizes.end(), OldXSizes.begin());
        std::copy(NewYSizes.begin(), NewYSizes.end(), OldYSizes.begin());
        std::copy(NewZSizes.begin(), NewZSizes.end(), OldZSizes.begin());
      }

    template<class ThreeDModelType>
    void CachedGravMagCalculator<ThreeDModelType>::CopyMeasPos(
        const MeasPosType &NewMeasPosX, const MeasPosType &NewMeasPosY,
        const MeasPosType &NewMeasPosZ)
      {
        OldMeasPosX.resize(NewMeasPosX.size());
        OldMeasPosY.resize(NewMeasPosY.size());
        OldMeasPosZ.resize(NewMeasPosZ.size());
        std::copy(NewMeasPosX.begin(), NewMeasPosX.end(), OldMeasPosX.begin());
        std::copy(NewMeasPosY.begin(), NewMeasPosY.end(), OldMeasPosY.begin());
        std::copy(NewMeasPosZ.begin(), NewMeasPosZ.end(), OldMeasPosZ.begin());
      }

    template<class ThreeDModelType>
    bool CachedGravMagCalculator<ThreeDModelType>::CheckMeasPosChange(
        const ThreeDModelType &Model)
      {
        bool change = true;
        // if all the sizes are the same as before then nothing has changed
        change = !(OldMeasPosX.size() == Model.GetMeasPosX().size()
            && OldMeasPosY.size() == Model.GetMeasPosY().size()
            && OldMeasPosZ.size() == Model.GetMeasPosZ().size());
        // if we already know that something has changed we do not need to perform the more expensive tests
        if (change)
          {
            //copy the new information into the cache
            CopyMeasPos(Model.GetMeasPosX(), Model.GetMeasPosY(), Model.GetMeasPosZ());
            return change;
          }
        //check whether any of the cell coordinates have changed
        bool xsame = std::equal(OldMeasPosX.begin(), OldMeasPosX.end(),
            Model.GetMeasPosX().begin());
        bool ysame = std::equal(OldMeasPosY.begin(), OldMeasPosY.end(),
            Model.GetMeasPosY().begin());
        bool zsame = std::equal(OldMeasPosZ.begin(), OldMeasPosZ.end(),
            Model.GetMeasPosZ().begin());
        //only if they are all the same we know that nothing has changed
        change = !(xsame && ysame && zsame);
        if (change)
          {
            CopyMeasPos(Model.GetMeasPosX(), Model.GetMeasPosY(), Model.GetMeasPosZ());
          }
        return change;
      }

    template<class ThreeDModelType>
    bool CachedGravMagCalculator<ThreeDModelType>::CheckGeometryChange(
        const ThreeDModelType &Model)
      {
        // by default we assume a change
        bool change = true;
        // if all the sizes are the same as before then nothing has changed
        change = (OldXSizes.size() != Model.GetXCellSizes().size()
            || OldYSizes.size() != Model.GetYCellSizes().size()
            || OldZSizes.size() != Model.GetZCellSizes().size());
        // if we already know that something has changed we do not need to perform the more expensive tests
        if (change)
          {
            //copy the new information into the cache
            CopySizes(Model.GetXCellSizes(), Model.GetYCellSizes(),
                Model.GetZCellSizes());
            return change;
          }
        //check whether any of the cell coordinates have changed
        bool xsame = std::equal(OldXSizes.begin(), OldXSizes.end(),
            Model.GetXCellSizes().begin());
        bool ysame = std::equal(OldYSizes.begin(), OldYSizes.end(),
            Model.GetYCellSizes().begin());
        bool zsame = std::equal(OldZSizes.begin(), OldZSizes.end(),
            Model.GetZCellSizes().begin());
        //only if they are all the same we know that nothing has changed
        change = !(xsame && ysame && zsame);
        if (change)
          {
            CopySizes(Model.GetXCellSizes(), Model.GetYCellSizes(),
                Model.GetZCellSizes());
          }
        return change;
      }

    template<class ThreeDModelType>
    bool CachedGravMagCalculator<ThreeDModelType>::CheckBackgroundChange(
        const ThreeDModelType &Model)
      {
        // by default we assume a change
        bool change = true;
        //check if either the size of the background densities or the thicknesses changed
        change = (OldBackgroundDens.size() != Model.GetBackgroundDensities().size()
            || OldBackgroundThick.size() != Model.GetBackgroundThicknesses().size());
        //check if one size changed copy the new values
        if (change)
          {
            OldBackgroundDens.resize(Model.GetBackgroundDensities().size(), 0.0);
            OldBackgroundThick.resize(Model.GetBackgroundThicknesses().size(), 0.0);
            std::copy(Model.GetBackgroundDensities().begin(),
                Model.GetBackgroundDensities().end(), OldBackgroundDens.begin());
            std::copy(Model.GetBackgroundThicknesses().begin(),
                Model.GetBackgroundThicknesses().end(), OldBackgroundThick.begin());
            return change;
          }
        //if the sizes are the same, we check whether the vectors still conatin the same values
        //for densities
        bool denssame = std::equal(OldBackgroundDens.begin(), OldBackgroundDens.end(),
            Model.GetBackgroundDensities().begin());
        //and for thickness
        bool thicksame = std::equal(OldBackgroundThick.begin(), OldBackgroundThick.end(),
            Model.GetBackgroundThicknesses().begin());
        change = !(denssame && thicksame);
        //if the content changed we copy the new values
        if (change)
          {
            OldBackgroundDens.resize(Model.GetBackgroundDensities().size(), 0.0);
            OldBackgroundThick.resize(Model.GetBackgroundThicknesses().size(), 0.0);
            std::copy(Model.GetBackgroundDensities().begin(),
                Model.GetBackgroundDensities().end(), OldBackgroundDens.begin());
            std::copy(Model.GetBackgroundThicknesses().begin(),
                Model.GetBackgroundThicknesses().end(), OldBackgroundThick.begin());
          }
        return change;
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
        //we have to do them separately to ensure that all functions are executed
        bool HasGeometryChanged = CheckGeometryChange(Model);
        bool HasMeasPosChanged = CheckMeasPosChange(Model);
        bool HasBGCHanged = CheckBackgroundChange(Model);
        bool HasChanged = HasGeometryChanged || HasMeasPosChanged || HasBGCHanged;
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
        bool HasChanged = CheckGeometryChange(Model) && CheckMeasPosChange(Model)
            && CheckBackgroundChange(Model);
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
