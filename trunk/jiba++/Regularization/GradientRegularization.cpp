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
                        XOperatorMatrix(index, ModelGeometry.IndexToOffset(i
                            + 1, j, k)) = 1.0;
                        XOperatorMatrix(index, index) = CenterValue;
                      }
                  }
              }
          }//end of for loop for x

        for (size_t i = 0; i < xsize; ++i)
          {
            for (size_t j = 0; j < ysize - 1; ++j)
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
                        YOperatorMatrix(index, ModelGeometry.IndexToOffset(i, j
                            + 1, k)) = 1.0;
                        YOperatorMatrix(index, index) = CenterValue;
                      }
                  }
              }
          }//end of for loop for x

        for (size_t i = 0; i < xsize; ++i)
          {
            for (size_t j = 0; j < ysize; ++j)
              {
                for (size_t k = 0; k < zsize - 1; ++k)
                  {
                    //the first index for the matrix element is always the same
                    const size_t index = ModelGeometry.IndexToOffset(i, j, k);
                    //we use forward differences for the gradient
                    //so we have two elements per matrix row
                    //for each operator matrix we have to check
                    //whether we want a tear for the cell boundary in that direction
                    if (TearModelZ.GetData()[i][j][k])
                      {
                        ZOperatorMatrix(index, ModelGeometry.IndexToOffset(i,
                            j, k + 1)) = 1.0;
                        ZOperatorMatrix(index, index) = CenterValue;
                      }
                  }
              }
          }//end of for loop for x
      }
  }
