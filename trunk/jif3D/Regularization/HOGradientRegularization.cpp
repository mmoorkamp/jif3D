//============================================================================
// Name        : GradientRegularization.cpp
// Author      : Jun 5, 2009
// Version     :
// Copyright   : 2009, mmoorkamp
//============================================================================

#include "HOGradientRegularization.h"
#include <cassert>
namespace jif3D
  {

    void HOGradientRegularization::ConstructOperator(
        const jif3D::ThreeDModelBase &ModelGeometry,
        const jif3D::ThreeDModelBase &TearModelX,
        const jif3D::ThreeDModelBase &TearModelY,
        const jif3D::ThreeDModelBase &TearModelZ)
      {
        const size_t xsize = ModelGeometry.GetModelShape()[0];
        const size_t ysize = ModelGeometry.GetModelShape()[1];
        const size_t zsize = ModelGeometry.GetModelShape()[2];

        const double a = 1.0 / 12.0;
        const double b = 2.0 / 3.0;
        //the inner part of the model where we can do the derivative in all directions
        //assign values to all three directions
        for (size_t i = 2; i < xsize - 2; ++i)
          {
            for (size_t j = 0; j < ysize; ++j)
              {
                for (size_t k = 0; k < zsize; ++k)
                  {
                    //the first index for the matrix element is always the same
                    const size_t index = ModelGeometry.IndexToOffset(i, j, k);
                    //we use forward differences for the gradient
                    //so we have two elements per matrix row
                    //for each operator matrix we have to check
                    //whether we want a tear for the cell boundary in that direction
                    if (TearModelX.GetData()[i][j][k])
                      {
                        XOperatorMatrix(index, ModelGeometry.IndexToOffset(i - 2, j, k)) =
                            a;
                        XOperatorMatrix(index, ModelGeometry.IndexToOffset(i - 1, j, k)) =
                            -b;
                        XOperatorMatrix(index, index) = Eps;
                        XOperatorMatrix(index, ModelGeometry.IndexToOffset(i + 2, j, k)) =
                            -a;
                        XOperatorMatrix(index, ModelGeometry.IndexToOffset(i + 1, j, k)) =
                            b;
                      }
                  }
              }
          } //end of for loop for x

        for (size_t i = 0; i < xsize; ++i)
          {
            for (size_t j = 2; j < ysize - 2; ++j)
              {
                for (size_t k = 0; k < zsize; ++k)
                  {
                    //the first index for the matrix element is always the same
                    const size_t index = ModelGeometry.IndexToOffset(i, j, k);
                    //we use forward differences for the gradient
                    //so we have two elements per matrix row
                    //for each operator matrix we have to check
                    //whether we want a tear for the cell boundary in that direction

                    if (TearModelY.GetData()[i][j][k])
                      {
                        YOperatorMatrix(index, ModelGeometry.IndexToOffset(i, j - 2, k)) =
                            a;
                        YOperatorMatrix(index, ModelGeometry.IndexToOffset(i, j - 1, k)) =
                            -b;
                        YOperatorMatrix(index, index) = Eps;
                        YOperatorMatrix(index, ModelGeometry.IndexToOffset(i, j + 2, k)) =
                            -a;
                        YOperatorMatrix(index, ModelGeometry.IndexToOffset(i, j + 1, k)) =
                            b;
                      }
                  }
              }
          } //end of for loop for y

        for (size_t i = 0; i < xsize; ++i)
          {
            for (size_t j = 0; j < ysize; ++j)
              {
                for (size_t k = 2; k < zsize - 2; ++k)
                  {
                    //the first index for the matrix element is always the same
                    const size_t index = ModelGeometry.IndexToOffset(i, j, k);
                    //we use forward differences for the gradient
                    //so we have two elements per matrix row
                    //for each operator matrix we have to check
                    //whether we want a tear for the cell boundary in that direction
                    if (TearModelZ.GetData()[i][j][k])
                      {
                        ZOperatorMatrix(index, ModelGeometry.IndexToOffset(i, j, k - 2)) =
                            a;
                        ZOperatorMatrix(index, ModelGeometry.IndexToOffset(i, j, k - 1)) =
                            -b;
                        ZOperatorMatrix(index, index) = Eps;
                        ZOperatorMatrix(index, ModelGeometry.IndexToOffset(i, j, k + 2)) =
                            -a;
                        ZOperatorMatrix(index, ModelGeometry.IndexToOffset(i, j, k + 1)) =
                            b;

                      }
                  }
              }
          } //end of for loop for x
      }
  }
