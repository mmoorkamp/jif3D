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
        const jiba::ThreeDModelBase &ModelGeometry)
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
                    const size_t index = ModelGeometry.IndexToOffset(i, j, k);
                    XOperatorMatrix(index, ModelGeometry.IndexToOffset(i + 1,
                        j, k)) = RightValue;
                    XOperatorMatrix(index, index) = CenterValue;
                    XOperatorMatrix(index, ModelGeometry.IndexToOffset(i - 1,
                        j, k)) = LeftValue;

                    YOperatorMatrix(index, ModelGeometry.IndexToOffset(i,
                        j + 1, k)) = RightValue;
                    YOperatorMatrix(index, index) = CenterValue;
                    YOperatorMatrix(index, ModelGeometry.IndexToOffset(i,
                        j - 1, k)) = LeftValue;

                    ZOperatorMatrix(index, ModelGeometry.IndexToOffset(i, j, k
                        + 1)) = RightValue;
                    ZOperatorMatrix(index, index) = CenterValue;
                    ZOperatorMatrix(index, ModelGeometry.IndexToOffset(i, j, k
                        - 1)) = LeftValue;
                  }
                //we can handle the border in z-direction within this loop
                //first the maximum z index
                size_t index = ModelGeometry.IndexToOffset(i, j, zsize - 1);
                XOperatorMatrix(index, ModelGeometry.IndexToOffset(i + 1, j,
                    zsize - 1)) = RightValue;
                XOperatorMatrix(index, index) = CenterValue;
                XOperatorMatrix(index, ModelGeometry.IndexToOffset(i - 1, j,
                    zsize - 1)) = LeftValue;

                YOperatorMatrix(index, ModelGeometry.IndexToOffset(i, j + 1,
                    zsize - 1)) = RightValue;
                YOperatorMatrix(index, index) = CenterValue;
                YOperatorMatrix(index, ModelGeometry.IndexToOffset(i, j - 1,
                    zsize - 1)) = LeftValue;

                //then the first z-index
                index = ModelGeometry.IndexToOffset(i, j, 0);
                XOperatorMatrix(index, ModelGeometry.IndexToOffset(i + 1, j, 0))
                    = RightValue;
                XOperatorMatrix(index, index) = CenterValue;
                XOperatorMatrix(index, ModelGeometry.IndexToOffset(i - 1, j, 0))
                    = LeftValue;

                YOperatorMatrix(index, ModelGeometry.IndexToOffset(i, j + 1, 0))
                    = RightValue;
                YOperatorMatrix(index, index) = CenterValue;
                YOperatorMatrix(index, ModelGeometry.IndexToOffset(i, j - 1, 0))
                    = LeftValue;
              }
          }
      }
  }
