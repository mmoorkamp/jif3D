//============================================================================
// Name        : ModelRefiner.cpp
// Author      : Feb 20, 2009
// Version     :
// Copyright   : 2009, mmoorkamp
//============================================================================

#include <numeric>
#include "ModelRefiner.h"
#include "../Global/FatalException.h"
namespace jiba
  {

    void ModelRefiner::RefineOneAxis(
        const ThreeDModelBase::t3DModelDim &OldCoordinates,
        const ThreeDModelBase::t3DModelDim &RefCoordinates,
        const ThreeDModelBase::t3DModelDim &OldSizes,
        ThreeDModelBase::t3DModelDim &NewSizes)
      {
        size_t oldindex = 0;
        size_t refineindex = 0;
        const size_t nold = OldCoordinates.size();
        const size_t nref = RefCoordinates.size();
        std::vector<double> TempCoord;
        for (size_t i = 0; i < nold; ++i)
          {
            TempCoord.push_back(OldCoordinates[i]);
            if (refineindex < nref && RefCoordinates[refineindex]
                == OldCoordinates[i])
              {
                ++refineindex;
              }
            const double CurrLimit = OldCoordinates[i] + OldSizes[i];

            while (refineindex < nref && RefCoordinates[refineindex]
                < CurrLimit)
              {
                TempCoord.push_back(RefCoordinates[refineindex]);
                ++refineindex;
              }

          }
        NewSizes.resize(boost::extents[TempCoord.size()]);
        std::adjacent_difference(TempCoord.begin(), TempCoord.end(),
            NewSizes.origin());
        std::rotate(NewSizes.begin(), NewSizes.begin() + 1, NewSizes.end());

        double lastsize = OldSizes[OldSizes.size() - 1];
        if (refineindex > 0)
          {
            lastsize += OldCoordinates[OldCoordinates.size() - 1] - TempCoord.back();
          }
        if (lastsize > 0.0)
          {
            NewSizes[NewSizes.size() - 1] = lastsize;
          }
        else
          {
            throw jiba::FatalException(
                "Negative cell size in model refinement !");
          }
      }

    void ModelRefiner::RefineAxes(const ThreeDModelBase &InputModel,
        ThreeDModelBase &RefinedModel)
      {

        RefineOneAxis(InputModel.GetXCoordinates(), RefiningXCoordinates,
            InputModel.GetXCellSizes(), RefinedModel.SetXCellSizes());
        RefineOneAxis(InputModel.GetYCoordinates(), RefiningYCoordinates,
            InputModel.GetYCellSizes(), RefinedModel.SetYCellSizes());
        RefineOneAxis(InputModel.GetZCoordinates(), RefiningZCoordinates,
            InputModel.GetZCellSizes(), RefinedModel.SetZCellSizes());
        RefinedModel .SetData().resize(
            boost::extents[RefinedModel.GetXCellSizes().size()][RefinedModel.GetYCellSizes().size()][RefinedModel.GetZCellSizes().size()]);
      }

    void FindEnd(size_t startindex, size_t &endindex, const double CellSize,
        const ThreeDModelBase::t3DModelDim &NewSizes)
      {
        double currsize = 0.0;
        endindex = startindex;
        while (currsize < CellSize && endindex < NewSizes.size())
          {
            currsize += NewSizes[endindex];
            ++endindex;
          }
      }

    void AssignValue(
        ThreeDModelBase::t3DModelData::array_view<3>::type &myview,
        const double value)
      {
        for (size_t i = 0; i < myview.shape()[0]; ++i)
          for (size_t j = 0; j < myview.shape()[1]; ++j)
            for (size_t k = 0; k < myview.shape()[2]; ++k)
              myview[i][j][k] = value;
      }

    void ModelRefiner::ProjectValues(const ThreeDModelBase &InputModel,
        ThreeDModelBase &RefinedModel)
      {
        const size_t oldxsize = InputModel.GetXCellSizes().size();
        const size_t oldysize = InputModel.GetYCellSizes().size();
        const size_t oldzsize = InputModel.GetZCellSizes().size();
        size_t startx = 0, endx = 0, starty = 0, endy = 0, startz = 0, endz = 0;
        for (size_t i = 0; i < oldxsize; ++i)
          {
            FindEnd(startx, endx, InputModel.GetXCellSizes()[i],
                RefinedModel.GetXCellSizes());
            for (size_t j = 0; j < oldysize; ++j)
              {
                FindEnd(starty, endy, InputModel.GetYCellSizes()[i],
                    RefinedModel.GetYCellSizes());
                for (size_t k = 0; k < oldzsize; ++k)
                  {
                    FindEnd(startz, endz, InputModel.GetZCellSizes()[i],
                        RefinedModel.GetZCellSizes());

                    typedef boost::multi_array_types::index_range range;
                    ThreeDModelBase::t3DModelData::array_view<3>::type myview =
                        RefinedModel.SetData()[boost::indices[range(startx,
                            endx)][range(starty, endy)][range(startz, endz)]];
                    AssignValue(myview, InputModel.GetData()[i][j][k]);
                    startz = endz;
                  }
                starty = endy;
                startz = 0;
              }
            startx = endx;
            starty = 0;
          }
      }

    ModelRefiner::ModelRefiner()
      {
        // TODO Auto-generated constructor stub

      }

    ModelRefiner::~ModelRefiner()
      {
        // TODO Auto-generated destructor stub
      }

  }
