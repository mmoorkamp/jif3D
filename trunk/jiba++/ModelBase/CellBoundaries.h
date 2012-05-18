//============================================================================
// Name        : CellBoundaries.h
// Author      : 18 May 2012
// Version     : 
// Copyright   : 2012, mm489
//============================================================================

#ifndef CELLBOUNDARIES_H_
#define CELLBOUNDARIES_H_

#include <algorithm>
#include "../Global/FatalException.h"
#include "../Global/convert.h"
#include "ThreeDModelBase.h"

namespace jiba
  {
    //! Given a coordinate within the modelling domain, return the index of the cell boundary closest to that coordinate.
    /*! Given a coordinate within the modelling domain, return the index of the cell boundary closest to that coordinate.
     * @param Coordinate The coordinate value, should lie within the modelling domain
     * @param CellBoundaries The array of cell boundary values, the first value is implicitly assumed to be zero
     * @return The index of the nearest boundary.
     */
    inline size_t FindNearestCellBoundary(double Coordinate,
        const jiba::ThreeDModelBase::t3DModelDim &CellBoundaries,
        const jiba::ThreeDModelBase::t3DModelDim &CellSizes)
      {
        size_t result = 0;
        const size_t ncells = CellBoundaries.num_elements();
        assert(CellBoundaries.size() == CellSizes.size());
        //we check whether the coordinate value is larger than the largest cell boundary
        //and throw an exception is true
        if (Coordinate > CellBoundaries[ncells - 1] + CellSizes[ncells - 1])
          throw jiba::FatalException(
              "Specified depth: " + jiba::stringify(Coordinate)
                  + " is outside domain limit: "
                  + jiba::stringify(CellBoundaries[ncells - 1]+ CellSizes[ncells - 1]));
        //we implictly assume the coordinates to start at zero
        //so the boundary values are the upper limits of the cells
        //and we search for the index of the next higher cell boundary
        while (Coordinate > (CellBoundaries[result] + CellSizes[result]) && result < ncells)
          {
            ++result;
          }
        //we do the substraction so that both values are positive, saves us an abs function call
        double uppdiff = Coordinate - CellBoundaries[result];
        double lowdiff = CellBoundaries[result] + CellSizes[result] - Coordinate;
        //return the index associated with the smaller positive difference
        if (uppdiff < lowdiff)
          return result;
        else
          return std::min(result + 1, ncells -1);
        //if we get here something is seriously wrong
        return -1;
      }
  }

#endif /* CELLBOUNDARIES_H_ */
