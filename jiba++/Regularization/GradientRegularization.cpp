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

    GradientRegularization::~GradientRegularization()
      {

      }

   void GradientRegularization::ConstructOperator(
        const jiba::ThreeDModelBase &ModelGeometry)
      {
        const size_t xsize = ModelGeometry.GetModelShape()[0];
        const size_t ysize = ModelGeometry.GetModelShape()[1];
        const size_t zsize = ModelGeometry.GetModelShape()[2];
        //the inner part of the matrix where we can do the derivative in all directions
        for (size_t i = 0; i < xsize - 1; ++i)
          {
            for (size_t j = 0; j < ysize - 1; ++j)
              {
                for (size_t k = 0; k < zsize - 1; ++k)
                  {
                    const size_t index = ModelGeometry.IndexToOffset(i, j, k);
                    XOperatorMatrix(index, ModelGeometry.IndexToOffset(i + 1,
                        j, k)) = 1.0;
                    XOperatorMatrix(index, index) = -1.0;

                    YOperatorMatrix(index, ModelGeometry.IndexToOffset(i,
                        j + 1, k)) = 1.0;
                    YOperatorMatrix(index, index) = -1.0;

                    ZOperatorMatrix(index, ModelGeometry.IndexToOffset(i, j, k
                        + 1)) = 1.0;
                    ZOperatorMatrix(index, index) = -1.0;
                  }
                //we can handle the border in z-direction within this loop
                const size_t index = ModelGeometry.IndexToOffset(i, j, zsize
                    - 1);
                XOperatorMatrix(index, ModelGeometry.IndexToOffset(i + 1, j,
                    zsize - 1)) = 1.0;
                XOperatorMatrix(index, index) = -1.0;

                YOperatorMatrix(index, ModelGeometry.IndexToOffset(i, j + 1,
                    zsize - 1)) = 1.0;
                YOperatorMatrix(index, index) = -1.0;
              }
          }

        //now take care of the border in x-direction
        for (size_t j = 0; j < ysize - 1; ++j)
          {
            for (size_t k = 0; k < zsize - 1; ++k)
              {
                const size_t index = ModelGeometry.IndexToOffset(xsize - 1, j,
                    k);
                YOperatorMatrix(index, ModelGeometry.IndexToOffset(xsize - 1, j
                    + 1, k)) = 1.0;
                YOperatorMatrix(index, index) = -1.0;
                ZOperatorMatrix(index, ModelGeometry.IndexToOffset(xsize - 1,
                    j, k + 1)) = 1.0;
                ZOperatorMatrix(index, index) = -1.0;
              }
          }

        //now take care of the border in y-direction
        for (size_t j = 0; j < xsize - 1; ++j)
          {
            for (size_t k = 0; k < zsize - 1; ++k)
              {
                const size_t index = ModelGeometry.IndexToOffset(j, ysize - 1,
                    k);
                XOperatorMatrix(index, ModelGeometry.IndexToOffset(j + 1, ysize
                    - 1, k)) = 1.0;
                XOperatorMatrix(index, index) = -1.0;
                ZOperatorMatrix(index, ModelGeometry.IndexToOffset(j,
                    ysize - 1, k + 1)) = 1.0;
                ZOperatorMatrix(index, index) = -1.0;
              }
          }
      }
  }
