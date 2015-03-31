//============================================================================
// Name        : CellBoundaries.cpp
// Author      : 4 Jul 2013
// Version     : 
// Copyright   : 2013, mm489
//============================================================================

#include <algorithm>
#include <cassert>
#include "../Global/FatalException.h"
#include "../Global/convert.h"
#include "CellBoundaries.h"

namespace jif3D
  {

    size_t FindNearestCellBoundary(double Coordinate,
        const jif3D::ThreeDModelBase::t3DModelDim &CellBoundaries,
        const jif3D::ThreeDModelBase::t3DModelDim &CellSizes)
      {
        size_t result = 0;
        const size_t ncells = CellBoundaries.num_elements();
        assert(CellBoundaries.size() == CellSizes.size());
        //we check whether the coordinate value is larger than the largest cell boundary
        //and throw an exception is true
        if (Coordinate > CellBoundaries[ncells - 1] + CellSizes[ncells - 1]
            || Coordinate < 0)
          throw jif3D::FatalException(
              "Specified depth: " + jif3D::stringify(Coordinate)
                  + " is outside domain limit: "
                  + jif3D::stringify(CellBoundaries[ncells - 1] + CellSizes[ncells - 1]),
              __FILE__, __LINE__);
        //we implictly assume the coordinates to start at zero
        //so the boundary values are the upper limits of the cells
        //and we search for the index of the next higher cell boundary
        while (Coordinate > (CellBoundaries[result] + CellSizes[result])
            && result < ncells)
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
          return std::min(result + 1, ncells - 1);
        //if we get here something is seriously wrong
        throw jif3D::FatalException("Unexpected error in FindNearestCellBoundary !",
            __FILE__, __LINE__);
      }

    size_t ConstructDepthIndices(std::vector<size_t> &MeasDepthIndices,
        std::vector<double> &ShiftDepth, const ThreeDModelBase &Model)
      {
        const size_t nmeas = Model.GetMeasPosX().size();
        size_t nlevels = 0;
        if (nmeas == 0)
          return nlevels;
        std::vector<size_t> ZIndices;

        for (size_t i = 0; i < nmeas; ++i)
          {
            size_t zindex = FindNearestCellBoundary(Model.GetMeasPosZ()[i],
                Model.GetZCoordinates(), Model.GetZCellSizes());
            std::vector<size_t>::iterator CurrIndex = std::find(ZIndices.begin(),
                ZIndices.end(), zindex);
            if (CurrIndex == ZIndices.end())
              {
                ZIndices.push_back(zindex);
                ShiftDepth.push_back(Model.GetZCoordinates()[zindex]);
                MeasDepthIndices.push_back(nlevels);
                ++nlevels;
              }
            else
              {
                MeasDepthIndices.push_back(std::distance(ZIndices.begin(), CurrIndex));
              }
          }
        return nlevels;
      }
  }

