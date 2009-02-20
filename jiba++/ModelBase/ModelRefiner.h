//============================================================================
// Name        : ModelRefiner.h
// Author      : Feb 20, 2009
// Version     :
// Copyright   : 2009, mmoorkamp
//============================================================================


#ifndef MODELREFINER_H_
#define MODELREFINER_H_
#include "ThreeDModelBase.h"
namespace jiba
  {

    class ModelRefiner
      {
    private:
      ThreeDModelBase::t3DModelDim RefiningXCoordinates;
      ThreeDModelBase::t3DModelDim RefiningYCoordinates;
      ThreeDModelBase::t3DModelDim RefiningZCoordinates;
      void RefineOneAxis(const ThreeDModelBase::t3DModelDim &OldCoordinates,
          const ThreeDModelBase::t3DModelDim &RefCoordinates,
          const ThreeDModelBase::t3DModelDim &OldSizes,
          ThreeDModelBase::t3DModelDim &NewSizes);
    public:
      void SetXCoordinates(ThreeDModelBase::t3DModelDim &XCoordinates)
        {
          RefiningXCoordinates.resize(boost::extents[XCoordinates.size()]);
          RefiningXCoordinates = XCoordinates;
        }
      void SetYCoordinates(ThreeDModelBase::t3DModelDim &YCoordinates)
        {
          RefiningYCoordinates.resize(boost::extents[YCoordinates.size()]);
          RefiningYCoordinates = YCoordinates;
        }
      void SetZCoordinates(ThreeDModelBase::t3DModelDim &ZCoordinates)
        {
          RefiningZCoordinates.resize(boost::extents[ZCoordinates.size()]);
          RefiningZCoordinates = ZCoordinates;
        }
      void RefineAxes(const ThreeDModelBase &InputModel,
          ThreeDModelBase &RefinedModel);
      void ProjectValues(const ThreeDModelBase &InputModel,
          ThreeDModelBase &RefinedModel);
      void RefineModel(const ThreeDModelBase &InputModel,
          ThreeDModelBase &RefinedModel)
        {
          RefineAxes(InputModel, RefinedModel);
          ProjectValues(InputModel, RefinedModel);
        }
      ModelRefiner();
      virtual ~ModelRefiner();
      };

  }

#endif /* MODELREFINER_H_ */
