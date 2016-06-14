//============================================================================
// Name        : GridOnlyModelCache.h
// Author      : 25 Oct 2013
// Version     : 
// Copyright   : 2013, mm489
//============================================================================

#ifndef GRIDONLYMODELCACHE_H_
#define GRIDONLYMODELCACHE_H_

namespace jif3D
  {
    template<class ThreeDModelType>
    class GridOnlyModelCache
      {
    public:
      typedef typename ThreeDModelType::t3DModelDim ModelDimType;
      typedef typename ThreeDModelType::tMeasPosVec MeasPosType;
    private:
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
    public:
      bool HasChanged(const ThreeDModelType &Model);
      GridOnlyModelCache();
      virtual ~GridOnlyModelCache();
      };

    template<class ThreeDModelType>
    GridOnlyModelCache<ThreeDModelType>::GridOnlyModelCache()
      {
        // TODO Auto-generated constructor stub

      }

    template<class ThreeDModelType>
    GridOnlyModelCache<ThreeDModelType>::~GridOnlyModelCache()
      {
        // TODO Auto-generated destructor stub
      }

    template<class ThreeDModelType>
    bool GridOnlyModelCache<ThreeDModelType>::HasChanged(const ThreeDModelType &Model)
      {
        bool HasGeometryChanged = CheckGeometryChange(Model);
        bool HasMeasPosChanged = CheckMeasPosChange(Model);
        bool HasChanged = HasGeometryChanged || HasMeasPosChanged;
        return HasChanged;
      }

    template<class ThreeDModelType>
    void GridOnlyModelCache<ThreeDModelType>::CopySizes(const ModelDimType &NewXSizes,
        const ModelDimType &NewYSizes, const ModelDimType &NewZSizes)
      {
        //make sure we have enough memory
        OldXSizes.resize(NewXSizes.size());
        OldYSizes.resize(NewYSizes.size());
        OldZSizes.resize(NewZSizes.size());
        //copy old sizes into member variables for next comparison
        std::copy(NewXSizes.begin(), NewXSizes.end(), OldXSizes.begin());
        std::copy(NewYSizes.begin(), NewYSizes.end(), OldYSizes.begin());
        std::copy(NewZSizes.begin(), NewZSizes.end(), OldZSizes.begin());
      }

    template<class ThreeDModelType>
    void GridOnlyModelCache<ThreeDModelType>::CopyMeasPos(const MeasPosType &NewMeasPosX,
        const MeasPosType &NewMeasPosY, const MeasPosType &NewMeasPosZ)
      {
        OldMeasPosX.resize(NewMeasPosX.size());
        OldMeasPosY.resize(NewMeasPosY.size());
        OldMeasPosZ.resize(NewMeasPosZ.size());
        std::copy(NewMeasPosX.begin(), NewMeasPosX.end(), OldMeasPosX.begin());
        std::copy(NewMeasPosY.begin(), NewMeasPosY.end(), OldMeasPosY.begin());
        std::copy(NewMeasPosZ.begin(), NewMeasPosZ.end(), OldMeasPosZ.begin());
      }

    template<class ThreeDModelType>
    bool GridOnlyModelCache<ThreeDModelType>::CheckMeasPosChange(
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
    bool GridOnlyModelCache<ThreeDModelType>::CheckGeometryChange(
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

  } /* namespace jif3D */
#endif /* GRIDONLYMODELCACHE_H_ */
