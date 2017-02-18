//============================================================================
// Name        : GradientRegularization.cpp
// Author      : Jun 5, 2009
// Version     :
// Copyright   : 2009, mmoorkamp
//============================================================================

#include "GradientRegularization.h"
#include <cassert>
namespace jif3D
  {

    void GradientRegularization::ConstructOperator(
        const jif3D::ThreeDModelBase &ModelGeometry,
        const jif3D::ThreeDModelBase &TearModelX,
        const jif3D::ThreeDModelBase &TearModelY,
        const jif3D::ThreeDModelBase &TearModelZ)
      {
        const size_t xsize = ModelGeometry.GetModelShape()[0];
        const size_t ysize = ModelGeometry.GetModelShape()[1];
        const size_t zsize = ModelGeometry.GetModelShape()[2];
        //we add a small number to the diagonal of the matrix
        //to make it non-singular
        const double CenterValue = -1.0 + Eps;
        const double AvgX = std::accumulate(ModelGeometry.GetXCellSizes().begin(),
            ModelGeometry.GetXCellSizes().end(), 0.0)
            / boost::numeric_cast<double>(ModelGeometry.GetXCellSizes().size());
        const double AvgY = std::accumulate(ModelGeometry.GetYCellSizes().begin(),
            ModelGeometry.GetYCellSizes().end(), 0.0)
            / boost::numeric_cast<double>(ModelGeometry.GetYCellSizes().size());
        const double AvgZ = std::accumulate(ModelGeometry.GetZCellSizes().begin(),
            ModelGeometry.GetZCellSizes().end(), 0.0)
            / boost::numeric_cast<double>(ModelGeometry.GetZCellSizes().size());
        //the inner part of the model where we can do the derivative in all directions
        //assign values to all three directions
        for (size_t i = 0; i < xsize - 1; ++i)
          {
            const double XCSize =
                WithCellsize ?
                    (ModelGeometry.GetXCellSizes()[i + 1]
                        + ModelGeometry.GetXCellSizes()[i]) / (2.0 * AvgX) :
                    1.0;
            for (size_t j = 0; j < ysize; ++j)
              {
                for (size_t k = 0; k < zsize; ++k)
                  {
                    //the first index for the matrix element is always the same
                    const size_t index = ModelGeometry.IndexToOffset(i, j, k);
                    //we use forward differences for the gradient
                    //so we have two elements per matrix row
                    //we use the tear model as a weight for the regularization
                    //so 1.0 means normal regularization, 0.0 no regularization
                    if (TearModelX.GetData()[i][j][k])
                      {
                        XOperatorMatrix(index, ModelGeometry.IndexToOffset(i + 1, j, k)) =
                            std::abs(TearModelX.GetData()[i][j][k]) / XCSize;
                        XOperatorMatrix(index, index) = CenterValue
                            * std::abs(TearModelX.GetData()[i][j][k]) / XCSize;
                      }
                  }
              }
          } //end of for loop for x

        for (size_t i = 0; i < xsize; ++i)
          {
            for (size_t j = 0; j < ysize - 1; ++j)
              {
                const double YCSize =
                    WithCellsize ?
                        (ModelGeometry.GetYCellSizes()[j + 1]
                            + ModelGeometry.GetYCellSizes()[j]) / (2.0 * AvgY) :
                        1.0;
                for (size_t k = 0; k < zsize; ++k)
                  {
                    //the first index for the matrix element is always the same
                    const size_t index = ModelGeometry.IndexToOffset(i, j, k);
                    //we use forward differences for the gradient
                    //so we have two elements per matrix row
                    //we use the tear model as a weight for the regularization
                    //so 1.0 means normal regularization, 0.0 no regularization
                    if (TearModelY.GetData()[i][j][k])
                      {
                        YOperatorMatrix(index, ModelGeometry.IndexToOffset(i, j + 1, k)) =
                            std::abs(TearModelY.GetData()[i][j][k]) / YCSize;
                        YOperatorMatrix(index, index) = CenterValue
                            * std::abs(TearModelY.GetData()[i][j][k]) / YCSize;
                      }
                  }
              } //end of for loop for x
          }
        for (size_t i = 0; i < xsize; ++i)
          {
            for (size_t j = 0; j < ysize; ++j)
              {
                for (size_t k = 0; k < zsize - 1; ++k)
                  {
                    const double ZCSize =
                        WithCellsize ?
                            (ModelGeometry.GetZCellSizes()[k + 1]
                                + ModelGeometry.GetZCellSizes()[k]) / (2.0 * AvgZ) :
                            1.0;
                    //the first index for the matrix element is always the same
                    const size_t index = ModelGeometry.IndexToOffset(i, j, k);
                    //we use forward differences for the gradient
                    //so we have two elements per matrix row
                    //we use the tear model as a weight for the regularization
                    //so 1.0 means normal regularization, 0.0 no regularization
                    if (TearModelZ.GetData()[i][j][k])
                      {
                        ZOperatorMatrix(index, ModelGeometry.IndexToOffset(i, j, k + 1)) =
                            std::abs(TearModelZ.GetData()[i][j][k]) / ZCSize;
                        ZOperatorMatrix(index, index) = CenterValue
                            * std::abs(TearModelZ.GetData()[i][j][k]) / ZCSize;
                      }
                  }
              }
          } //end of for loop for x
      }
  }
