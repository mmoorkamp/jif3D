//============================================================================
// Name        : GridBackgroundModelCache.h
// Author      : 25 Oct 2013
// Version     :
// Copyright   : 2013, mm489
//============================================================================

#ifndef GRIDBACKGROUNDMODELCACHE_H_
#define GRIDBACKGROUNDMODELCACHE_H_

#include "../Global/Jif3DGlobal.h"
#include "GridOnlyModelCache.h"

namespace jif3D
  {
    template<class ThreeDModelType>
    class J3DEXPORT GridBackgroundModelCache
      {
    public:
      typedef typename ThreeDModelType::tBackgroundVec BGType;
    private:
      BGType OldBackgroundDens;
      BGType OldBackgroundThick;
      GridOnlyModelCache<ThreeDModelType> GridCache;
      //! Check wether the 1D background has changed since the last calculation
      bool CheckBackgroundChange(const ThreeDModelType &Model);
    public:
      bool HasChanged(const ThreeDModelType &Model);
      GridBackgroundModelCache()
        {
        }
      virtual ~GridBackgroundModelCache()
        {
        }
      };

    template<class ThreeDModelType>
    bool GridBackgroundModelCache<ThreeDModelType>::HasChanged(
        const ThreeDModelType &Model)
      {
        bool GridChange = GridCache.HasChanged(Model);
        bool bgChange = CheckBackgroundChange(Model);
        return GridChange || bgChange;
      }

    template<class ThreeDModelType>
    bool GridBackgroundModelCache<ThreeDModelType>::CheckBackgroundChange(
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
        //if the sizes are the same, we check whether the vectors still contain the same values
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

  } /* namespace jif3D */
#endif /* GRIDBACKGROUNDMODELCACHE_H_ */
