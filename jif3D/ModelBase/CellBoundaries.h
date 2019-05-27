//============================================================================
// Name        : CellBoundaries.h
// Author      : 18 May 2012
// Version     :
// Copyright   : 2012, mm489
//============================================================================

#ifndef CELLBOUNDARIES_H_
#define CELLBOUNDARIES_H_

#include "../Global/Jif3DGlobal.h"
#include "ThreeDModelBase.h"

namespace jif3D
  {
    /** \addtogroup modelbase Basic classes and routines for 3D models */
    /* @{ */
    //! Given a coordinate within the modelling domain, return the index of the cell boundary closest to that coordinate.
    /*! Given a coordinate within the modelling domain, return the index of the cell boundary closest to that coordinate.
     * @param Coordinate The coordinate value, should lie within the modelling domain
     * @param CellBoundaries The array of cell boundary values, the first value is implicitly assumed to be zero
     * @return The index of the nearest boundary.
     */
  J3DEXPORT size_t FindNearestCellBoundary(double Coordinate,
        const jif3D::ThreeDModelBase::t3DModelDim &CellBoundaries,
        const jif3D::ThreeDModelBase::t3DModelDim &CellSizes);

  J3DEXPORT size_t ConstructDepthIndices(std::vector<size_t> &MeasDepthIndices,
        std::vector<double> &ShiftDepth, const ThreeDModelBase &Model, const std::vector<double> &MeasPosZ);
  /* @} */
  }

#endif /* CELLBOUNDARIES_H_ */
