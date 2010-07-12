//============================================================================
// Name        : GradientRegularization.cpp
// Author      : Jun 5, 2009
// Version     :
// Copyright   : 2009, mmoorkamp
//============================================================================


#include "GradientRegularization.h"
#include <cassert>
namespace jiba
  {


    void GradientRegularization::ConstructOperator(
        const jiba::ThreeDModelBase &ModelGeometry,
        const jiba::ThreeDModelBase &TearModelX,
        const jiba::ThreeDModelBase &TearModelY,
        const jiba::ThreeDModelBase &TearModelZ)
      {
        const size_t xsize = ModelGeometry.GetModelShape()[0];
        const size_t ysize = ModelGeometry.GetModelShape()[1];
        const size_t zsize = ModelGeometry.GetModelShape()[2];
        //we add a small number to the diagonal of the matrix
        //to make it non-singular
        const double CenterValue = -1.0 + Eps;
        //the inner part of the model where we can do the derivative in all directions
        //assign values to all three directions
        for (size_t i = 0; i < xsize - 1; ++i)
          {
            for (size_t j = 0; j < ysize - 1; ++j)
              {
                for (size_t k = 0; k < zsize - 1; ++k)
                  {
                    //the first index for the matrix element is always the same
                    const size_t index = ModelGeometry.IndexToOffset(i, j, k);
                    //we use forward differences for the gradient
                    //so we have two elements per matrix row
                    if (TearModelX.GetData()[i][j][k])
                      {
                        XOperatorMatrix(index, ModelGeometry.IndexToOffset(i
                            + 1, j, k)) = 1.0;
                        XOperatorMatrix(index, index) = CenterValue;
                      }
                    if (TearModelY.GetData()[i][j][k])
                      {
                        YOperatorMatrix(index, ModelGeometry.IndexToOffset(i, j
                            + 1, k)) = 1.0;
                        YOperatorMatrix(index, index) = CenterValue;
                      }
                    if (TearModelZ.GetData()[i][j][k])
                      {
                        ZOperatorMatrix(index, ModelGeometry.IndexToOffset(i,
                            j, k + 1)) = 1.0;
                        ZOperatorMatrix(index, index) = CenterValue;
                      }
                  }

                //we can handle the border in z-direction within this loop
                //there we only use the gradient in x-direction and y-direction
                const size_t index = ModelGeometry.IndexToOffset(i, j, zsize
                    - 1);
                if (TearModelX.GetData()[i][j][zsize - 1])
                  {
                    XOperatorMatrix(index, ModelGeometry.IndexToOffset(i + 1,
                        j, zsize - 1)) = 1.0;
                    XOperatorMatrix(index, index) = CenterValue;
                  }
                if (TearModelY.GetData()[i][j][zsize - 1])
                  {
                    YOperatorMatrix(index, ModelGeometry.IndexToOffset(i,
                        j + 1, zsize - 1)) = 1.0;
                    YOperatorMatrix(index, index) = CenterValue;
                  }
              }
          }

        //now take care of the border in x-direction
        //there we only use the gradient in y-direction and z-direction
        for (size_t j = 0; j < ysize - 1; ++j)
          {
            for (size_t k = 0; k < zsize - 1; ++k)
              {
                const size_t index = ModelGeometry.IndexToOffset(xsize - 1, j,
                    k);
                if (TearModelY.GetData()[xsize - 1][j][k])
                  {
                    YOperatorMatrix(index, ModelGeometry.IndexToOffset(xsize
                        - 1, j + 1, k)) = 1.0;
                    YOperatorMatrix(index, index) = CenterValue;
                  }
                if (TearModelZ.GetData()[xsize - 1][j][k])
                  {
                    ZOperatorMatrix(index, ModelGeometry.IndexToOffset(xsize
                        - 1, j, k + 1)) = 1.0;
                    ZOperatorMatrix(index, index) = CenterValue;
                  }
              }
          }

        //now take care of the border in y-direction
        //there we only use the gradient in x-direction and z-direction
        for (size_t j = 0; j < xsize - 1; ++j)
          {
            for (size_t k = 0; k < zsize - 1; ++k)
              {
                const size_t index = ModelGeometry.IndexToOffset(j, ysize - 1,
                    k);
                if (TearModelX.GetData()[j][ysize - 1][k])
                  {
                    XOperatorMatrix(index, ModelGeometry.IndexToOffset(j + 1,
                        ysize - 1, k)) = 1.0;
                    XOperatorMatrix(index, index) = CenterValue;
                  }
                if (TearModelZ.GetData()[j][ysize - 1][k])
                  {
                    ZOperatorMatrix(index, ModelGeometry.IndexToOffset(j, ysize
                        - 1, k + 1)) = 1.0;
                    ZOperatorMatrix(index, index) = CenterValue;
                  }
              }
          }
      }
  }
