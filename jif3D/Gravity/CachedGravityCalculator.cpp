//============================================================================
// Name        : CachedGravityCalculator.cpp
// Author      : Dec 12, 2008
// Version     :
// Copyright   : 2008, mmoorkamp
//============================================================================


#include "CachedGravityCalculator.h"

namespace jif3D
  {

    CachedGravityCalculator::CachedGravityCalculator(boost::shared_ptr<
        ThreeDGravityImplementation> TheImp) :
      ThreeDGravityCalculator(TheImp), HaveCache(false), OldXSizes(),
          OldYSizes(), OldZSizes(), OldMeasPosX(), OldMeasPosY(),
          OldMeasPosZ(), OldBackgroundThick(), OldBackgroundDens()
      {

      }

    CachedGravityCalculator::~CachedGravityCalculator()
      {
      }

    void CachedGravityCalculator::CopySizes(
        const ThreeDGravityModel::t3DModelDim &NewXSizes,
        const ThreeDGravityModel::t3DModelDim &NewYSizes,
        const ThreeDGravityModel::t3DModelDim &NewZSizes)
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

    void CachedGravityCalculator::CopyMeasPos(
        const ThreeDGravityModel::tMeasPosVec &NewMeasPosX,
        const ThreeDGravityModel::tMeasPosVec &NewMeasPosY,
        const ThreeDGravityModel::tMeasPosVec &NewMeasPosZ)
      {
        OldMeasPosX.resize(NewMeasPosX.size());
        OldMeasPosY.resize(NewMeasPosY.size());
        OldMeasPosZ.resize(NewMeasPosZ.size());
        std::copy(NewMeasPosX.begin(), NewMeasPosX.end(), OldMeasPosX.begin());
        std::copy(NewMeasPosY.begin(), NewMeasPosY.end(), OldMeasPosY.begin());
        std::copy(NewMeasPosZ.begin(), NewMeasPosZ.end(), OldMeasPosZ.begin());
      }

    bool CachedGravityCalculator::CheckMeasPosChange(
        const ThreeDGravityModel &Model)
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
            CopyMeasPos(Model.GetMeasPosX(), Model.GetMeasPosY(),
                Model.GetMeasPosZ());
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
            CopyMeasPos(Model.GetMeasPosX(), Model.GetMeasPosY(),
                Model.GetMeasPosZ());
          }
        return change;
      }

    bool CachedGravityCalculator::CheckGeometryChange(
        const ThreeDGravityModel &Model)
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

    bool CachedGravityCalculator::CheckBackgroundChange(
        const ThreeDGravityModel &Model)
      {
        // by default we assume a change
        bool change = true;
        //check if either the size of the background densities or the thicknesses changed
        change = (OldBackgroundDens.size()
            != Model.GetBackgroundDensities().size()
            || OldBackgroundThick.size()
                != Model.GetBackgroundThicknesses().size());
        //check if one size changed copy the new values
        if (change)
          {
            OldBackgroundDens.resize(Model.GetBackgroundDensities().size(), 0.0);
            OldBackgroundThick.resize(Model.GetBackgroundThicknesses().size(),
                0.0);
            std::copy(Model.GetBackgroundDensities().begin(),
                Model.GetBackgroundDensities().end(), OldBackgroundDens.begin());
            std::copy(Model.GetBackgroundThicknesses().begin(),
                Model.GetBackgroundThicknesses().end(),
                OldBackgroundThick.begin());
            return change;
          }
        //if the sizes are the same, we check whether the vectors still conatin the same values
        //for densities
        bool denssame = std::equal(OldBackgroundDens.begin(),
            OldBackgroundDens.end(), Model.GetBackgroundDensities().begin());
        //and for thickness
        bool thicksame = std::equal(OldBackgroundThick.begin(),
            OldBackgroundThick.end(), Model.GetBackgroundThicknesses().begin());
        change = !(denssame && thicksame);
        //if the content changed we copy the new values
        if (change)
          {
            OldBackgroundDens.resize(Model.GetBackgroundDensities().size(), 0.0);
            OldBackgroundThick.resize(Model.GetBackgroundThicknesses().size(),
                0.0);
            std::copy(Model.GetBackgroundDensities().begin(),
                Model.GetBackgroundDensities().end(), OldBackgroundDens.begin());
            std::copy(Model.GetBackgroundThicknesses().begin(),
                Model.GetBackgroundThicknesses().end(),
                OldBackgroundThick.begin());
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
    rvec CachedGravityCalculator::Calculate(const ThreeDGravityModel &Model)
      {
        //check that all modeling information is consistent
        CheckModelConsistency(Model);
        //check whether the model geometry has changed since the last calculation
        //we have to do them separately to ensure that all functions are executed
        bool HasGeometryChanged = CheckGeometryChange(Model);
        bool HasMeasPosChanged = CheckMeasPosChange(Model);
        bool HasBGCHanged = CheckBackgroundChange(Model);
        bool HasChanged = HasGeometryChanged || HasMeasPosChanged
            || HasBGCHanged;
        //if we have cached information and nothing has changed, we use the cache
        if (HaveCache && !HasChanged)
          {
            return CalculateCachedResult(Model);
          }
        //we only get here if we need to recalculate
        //we have to make sure the calculation finishes properly before we can guarantee the cache and return the result
        SetCurrentSensitivities().resize(Imp.get()->RawDataPerMeasurement(),
            Model.GetDensities().num_elements()
                + Model.GetBackgroundDensities().size());
        rvec result = CalculateNewModel(Model);
        HaveCache = true;

        return result;
      }

    rvec CachedGravityCalculator::LQDerivative(const ThreeDGravityModel &Model,
        const rvec &Misfit)
      {
        //check that all modeling information is consistent
        CheckModelConsistency(Model);
        //check whether the model geometry has changed since the last calculation
        bool HasChanged = CheckGeometryChange(Model) && CheckMeasPosChange(
            Model) && CheckBackgroundChange(Model);
        //if we have cached information and nothing has changed, we use the cache
        if (HaveCache && !HasChanged)
          {
            return CachedLQDerivative(Model, Misfit);
          }
        //we only get here if we need to recalculate
        //we have to make sure the calculation finishes properly before we can guarantee the cache and return the result
        SetCurrentSensitivities().resize(Imp.get()->GetDataPerMeasurement(),
            Model.GetDensities().num_elements()
                + Model.GetBackgroundDensities().size());
        return ThreeDGravityCalculator::LQDerivative(Model, Misfit);

      }
  }
