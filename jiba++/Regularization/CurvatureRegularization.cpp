//============================================================================
// Name        : CurvatureRegularization.cpp
// Author      : Jan 11, 2010
// Version     : 
// Copyright   : 2010, mmoorkamp
//============================================================================


#include "CurvatureRegularization.h"

namespace jiba
  {

    CurvatureRegularization::~CurvatureRegularization()
      {

      }

    void CurvatureRegularization::ConstructOperator(
        const jiba::ThreeDModelBase &ModelGeometry,
        const jiba::ThreeDModelBase &TearModelX,
        const jiba::ThreeDModelBase &TearModelY,
        const jiba::ThreeDModelBase &TearModelZ)
      {
        const size_t xsize = ModelGeometry.GetModelShape()[0];
        const size_t ysize = ModelGeometry.GetModelShape()[1];
        const size_t zsize = ModelGeometry.GetModelShape()[2];
        const double CenterValue = -2.0 + Eps;
        const double LeftValue = 1.0;
        const double RightValue = 1.0;
        //the inner part of the matrix where we can do the derivative in all directions
        for (size_t i = 1; i < xsize - 1; ++i)
          {
            for (size_t j = 1; j < ysize - 1; ++j)
              {
                for (size_t k = 1; k < zsize - 1; ++k)
                  {
                    //save the storage index of the current element in the inversion domain
                    const size_t index = ModelGeometry.IndexToOffset(i, j, k);
                    //for each direction we have to check whether we want to introduce a tear
                    //we have to check the current element and the previous one,
                    //as curvature based regularization depends on elements on both
                    //sides of the current element, if both are different from zero
                    //we add the appropriate values to the operator matrix
                    if (TearModelX.GetData()[i][j][k] && TearModelX.GetData()[i
                        - 1][j][k])
                      {
                        XOperatorMatrix(index, ModelGeometry.IndexToOffset(i
                            + 1, j, k)) = RightValue;
                        XOperatorMatrix(index, index) = CenterValue;
                        XOperatorMatrix(index, ModelGeometry.IndexToOffset(i
                            - 1, j, k)) = LeftValue;
                      }
                    if (TearModelY.GetData()[i][j][k]
                        && TearModelY.GetData()[i][j - 1][k])
                      {
                        YOperatorMatrix(index, ModelGeometry.IndexToOffset(i, j
                            + 1, k)) = RightValue;
                        YOperatorMatrix(index, index) = CenterValue;
                        YOperatorMatrix(index, ModelGeometry.IndexToOffset(i, j
                            - 1, k)) = LeftValue;
                      }
                    if (TearModelZ.GetData()[i][j][k]
                        && TearModelZ.GetData()[i][j][k - 1])
                      {
                        ZOperatorMatrix(index, ModelGeometry.IndexToOffset(i,
                            j, k + 1)) = RightValue;
                        ZOperatorMatrix(index, index) = CenterValue;
                        ZOperatorMatrix(index, ModelGeometry.IndexToOffset(i,
                            j, k - 1)) = LeftValue;
                      }
                  }
                //we can handle the border in z-direction within this loop
                //first the maximum z index
                size_t index = ModelGeometry.IndexToOffset(i, j, zsize - 1);
                if (TearModelX.GetData()[i][j][zsize - 1]
                    && TearModelX.GetData()[i - 1][j][zsize - 1])
                  {
                    XOperatorMatrix(index, ModelGeometry.IndexToOffset(i + 1,
                        j, zsize - 1)) = RightValue;
                    XOperatorMatrix(index, index) = CenterValue;
                    XOperatorMatrix(index, ModelGeometry.IndexToOffset(i - 1,
                        j, zsize - 1)) = LeftValue;
                  }
                if (TearModelY.GetData()[i][j][zsize - 1]
                    && TearModelY.GetData()[i][j - 1][zsize - 1])
                  {
                    YOperatorMatrix(index, ModelGeometry.IndexToOffset(i,
                        j + 1, zsize - 1)) = RightValue;
                    YOperatorMatrix(index, index) = CenterValue;
                    YOperatorMatrix(index, ModelGeometry.IndexToOffset(i,
                        j - 1, zsize - 1)) = LeftValue;
                  }
                //then the first z-index
                if (TearModelX.GetData()[i][j][0]
                    && TearModelX.GetData()[i - 1][j][0])
                  {
                    index = ModelGeometry.IndexToOffset(i, j, 0);
                    XOperatorMatrix(index, ModelGeometry.IndexToOffset(i + 1,
                        j, 0)) = RightValue;
                    XOperatorMatrix(index, index) = CenterValue;
                    XOperatorMatrix(index, ModelGeometry.IndexToOffset(i - 1,
                        j, 0)) = LeftValue;
                  }
                if (TearModelY.GetData()[i][j][0] && TearModelY.GetData()[i][j
                    - 1][0])
                  {
                    YOperatorMatrix(index, ModelGeometry.IndexToOffset(i,
                        j + 1, 0)) = RightValue;
                    YOperatorMatrix(index, index) = CenterValue;
                    YOperatorMatrix(index, ModelGeometry.IndexToOffset(i,
                        j - 1, 0)) = LeftValue;
                  }
              }
          }
      }
  }
