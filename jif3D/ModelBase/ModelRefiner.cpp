//============================================================================
// Name        : ModelRefiner.cpp
// Author      : Feb 20, 2009
// Version     :
// Copyright   : 2009, mmoorkamp
//============================================================================

#include <numeric>
#include "ModelRefiner.h"
#include "../Global/FatalException.h"

 #pragma GCC diagnostic ignored "-Wuninitialized"

namespace jif3D
  {

    void ModelRefiner::RefineOneAxis(
        const ThreeDModelBase::t3DModelDim &OldCoordinates,
        const ThreeDModelBase::t3DModelDim &RefCoordinates,
        const ThreeDModelBase::t3DModelDim &OldSizes,
        ThreeDModelBase::t3DModelDim &NewSizes)
      {
        //we need an index for the current position on the refined axes
        size_t refineindex = 0;
        const size_t nold = OldCoordinates.size();
        const size_t nref = RefCoordinates.size();
        //a vector is slightly easier to handle than boost::multi_array, so we work with it and copy the result
        std::vector<double> TempCoord;
        // we go through all the coordinates of the original axes
        for (size_t i = 0; i < nold; ++i)
          {
            //we guarantee that the original coordinates are also in the refined coordinates
            TempCoord.push_back(OldCoordinates[i]);
            //we do not want to have a coordinate twice
            //even if the original coordinates also appear in the refined coordinates
            if (refineindex < nref && RefCoordinates[refineindex]
                == OldCoordinates[i])
              {
                ++refineindex;
              }
            //this is the end of the current cell
            const double CurrLimit = OldCoordinates[i] + OldSizes[i];
            // as long as we haven't exceeded the index of the refinement coordinates
            //and we are still in the current cell
            while (refineindex < nref && RefCoordinates[refineindex]
                < CurrLimit)
              {
                //we add the refinement coordinates
                TempCoord.push_back(RefCoordinates[refineindex]);
                ++refineindex;
              }

          }
        //we work with sizes so we transform the coordinates
        //to cell sizes
        NewSizes.resize(boost::extents[TempCoord.size()]);
        std::adjacent_difference(TempCoord.begin(), TempCoord.end(),
            NewSizes.origin());
        std::rotate(NewSizes.begin(), NewSizes.begin() + 1, NewSizes.end());
        //finally adjust the size of the last cell so that the modeling domain
        //is the same as before
        double lastsize = OldSizes[OldSizes.size() - 1];
        if (refineindex > 0)
          {
            lastsize += OldCoordinates[OldCoordinates.size() - 1]
                - TempCoord.back();
          }
        if (lastsize > 0.0)
          {
            NewSizes[NewSizes.size() - 1] = lastsize;
          }
        else
          {
            throw jif3D::FatalException(
                "Negative cell size in model refinement !", __FILE__, __LINE__);
          }
      }

    void ModelRefiner::RefineAxes(const ThreeDModelBase &InputModel,
        ThreeDModelBase &RefinedModel)
      {
        //Refine all three axes and allocate memory for the new grid
        RefineOneAxis(InputModel.GetXCoordinates(), RefiningXCoordinates,
            InputModel.GetXCellSizes(), RefinedModel.SetXCellSizes());
        RefineOneAxis(InputModel.GetYCoordinates(), RefiningYCoordinates,
            InputModel.GetYCellSizes(), RefinedModel.SetYCellSizes());
        RefineOneAxis(InputModel.GetZCoordinates(), RefiningZCoordinates,
            InputModel.GetZCellSizes(), RefinedModel.SetZCellSizes());
        RefinedModel .SetData().resize(
            boost::extents[RefinedModel.GetXCellSizes().size()][RefinedModel.GetYCellSizes().size()][RefinedModel.GetZCellSizes().size()]);
      }
    //find the index of the last refined cell that still corresponds to the current old cell
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
    //assign the same value to the subset of refined cells that correspond to one coarse cell
    void AssignValue(
        ThreeDModelBase::t3DModelData::array_view<3>::type &myview,
        const double value)
      {
        for (size_t i = 0; i < myview.shape()[0]; ++i)
          for (size_t j = 0; j < myview.shape()[1]; ++j)
            for (size_t k = 0; k < myview.shape()[2]; ++k)
              myview[i][j][k] = value;
      }

    //assign the same value to the subset of refined cells that correspond to one coarse cell
    double CombineValues(
        ThreeDModelBase::t3DModelData::array_view<3>::type &myview)
      {
        double result = 0;
        for (size_t i = 0; i < myview.shape()[0]; ++i)
          for (size_t j = 0; j < myview.shape()[1]; ++j)
            for (size_t k = 0; k < myview.shape()[2]; ++k)
              result += myview[i][j][k];
        return result;
      }

    void ModelRefiner::ProjectValues(const ThreeDModelBase &InputModel,
        ThreeDModelBase &RefinedModel)
      {
        //save the original size of the grid
        const size_t oldxsize = InputModel.GetXCellSizes().size();
        const size_t oldysize = InputModel.GetYCellSizes().size();
        const size_t oldzsize = InputModel.GetZCellSizes().size();
        //we have to keep track of the current start und end index
        //of the original coarse grid cell in the new refined grid
        //in all three directions
        size_t startx = 0, endx = 0, starty = 0, endy = 0, startz = 0, endz = 0;
        //we loop through all old cells
        for (size_t i = 0; i < oldxsize; ++i)
          {
            //find the end of the current cell in x-direction in the refined grid
            FindEnd(startx, endx, InputModel.GetXCellSizes()[i],
                RefinedModel.GetXCellSizes());
            for (size_t j = 0; j < oldysize; ++j)
              {
                //find the end of the current cell in y-direction in the refined grid
                FindEnd(starty, endy, InputModel.GetYCellSizes()[j],
                    RefinedModel.GetYCellSizes());
                for (size_t k = 0; k < oldzsize; ++k)
                  {
                    //find the end of the current cell in z-direction in the refined grid
                    FindEnd(startz, endz, InputModel.GetZCellSizes()[k],
                        RefinedModel.GetZCellSizes());
                    //now we know the indices of the old cell in the new grid
                    //and we can assign the right value to these cells
                    typedef boost::multi_array_types::index_range range;
                    ThreeDModelBase::t3DModelData::array_view<3>::type myview =
                        RefinedModel.SetData()[boost::indices[range(startx,
                            endx)][range(starty, endy)][range(startz, endz)]];
                    AssignValue(myview, InputModel.GetData()[i][j][k]);
                    //the next cell starts at the end of the current cell
                    startz = endz;
                  }
                starty = endy;
                //we start at the top again
                startz = 0;
              }
            startx = endx;
            starty = 0;
          }
      }

    jif3D::rvec ModelRefiner::CombineGradient(const jif3D::rvec &FineGradient,
        const ThreeDModelBase &CoarseModel, const ThreeDModelBase &RefinedModel)
      {
        //save the original size of the grid
        const size_t oldxsize = CoarseModel.GetXCellSizes().size();
        const size_t oldysize = CoarseModel.GetYCellSizes().size();
        const size_t oldzsize = CoarseModel.GetZCellSizes().size();
        //we have to keep track of the current start und end index
        //of the original coarse grid cell in the new refined grid
        //in all three directions
        ThreeDModelBase GradientModel(RefinedModel);
        jif3D::rvec CoarseGradient(CoarseModel.GetData().num_elements());
        std::copy(FineGradient.begin(), FineGradient.end(),
            GradientModel.SetData().origin());
        size_t startx = 0, endx = 0, starty = 0, endy = 0, startz = 0, endz = 0;
        //we loop through all old cells
        for (size_t i = 0; i < oldxsize; ++i)
          {
            //find the end of the current cell in x-direction in the refined grid
            FindEnd(startx, endx, CoarseModel.GetXCellSizes()[i],
                RefinedModel.GetXCellSizes());
            for (size_t j = 0; j < oldysize; ++j)
              {
                //find the end of the current cell in y-direction in the refined grid
                FindEnd(starty, endy, CoarseModel.GetYCellSizes()[j],
                    RefinedModel.GetYCellSizes());
                for (size_t k = 0; k < oldzsize; ++k)
                  {
                    //find the end of the current cell in z-direction in the refined grid
                    FindEnd(startz, endz, CoarseModel.GetZCellSizes()[k],
                        RefinedModel.GetZCellSizes());
                    typedef boost::multi_array_types::index_range range;
                    ThreeDModelBase::t3DModelData::array_view<3>::type myview =
                        GradientModel.SetData()[boost::indices[range(startx,
                            endx)][range(starty, endy)][range(startz, endz)]];
                    CoarseGradient(CoarseModel.IndexToOffset(i, j, k))
                        = CombineValues(myview);
                    //the next cell starts at the end of the current cell
                    startz = endz;
                  }
                starty = endy;
                //we start at the top again
                startz = 0;
              }
            startx = endx;
            starty = 0;
          }
        return CoarseGradient;
      }

    ModelRefiner::ModelRefiner() :
      RefiningXCoordinates(), RefiningYCoordinates(), RefiningZCoordinates()
      {

      }

    ModelRefiner::~ModelRefiner()
      {

      }

  }
