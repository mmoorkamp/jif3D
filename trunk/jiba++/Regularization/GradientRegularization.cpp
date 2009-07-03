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

    GradientRegularization::comp_mat GradientRegularization::SmoothnessOperator(
        const jiba::ThreeDModelBase &ModelGeometry)
      {
        comp_mat Result(ModelGeometry.GetData().num_elements(),
            ModelGeometry.GetData().num_elements());
        const size_t xsize = ModelGeometry.GetData().shape()[0];
        const size_t ysize = ModelGeometry.GetData().shape()[1];
        const size_t zsize = ModelGeometry.GetData().shape()[2];
        //the inner part of the matrix where we can do the derivative in all directions
        for (size_t i = 0; i < xsize - 1; ++i)
          {
            for (size_t j = 0; j < ysize - 1; ++j)
              {
                for (size_t k = 0; k < zsize - 1; ++k)
                  {
                    const size_t index = ModelGeometry.IndexToOffset(i, j, k);
                    Result(index, ModelGeometry.IndexToOffset(i + 1, j, k))
                        = 1.0;
                    Result(index, index) = -3.0;
                    Result(index, ModelGeometry.IndexToOffset(i, j + 1, k))
                        = 1.0;
                    Result(index, ModelGeometry.IndexToOffset(i, j, k + 1))
                        = 1.0;
                  }
              }
          }
        //now take care of the border in x-direction
        for (size_t j = 0; j < ysize - 1; ++j)
          {
            for (size_t k = 0; k < zsize - 1; ++k)
              {
                const size_t index = ModelGeometry.IndexToOffset(xsize - 1, j,
                    k);
                Result(index, index) = -2.0;
                Result(index, ModelGeometry.IndexToOffset(xsize - 1, j + 1, k))
                    = 1.0;
                Result(index, ModelGeometry.IndexToOffset(xsize - 1, j, k + 1))
                    = 1.0;
              }
          }

        //now take care of the border in y-direction
        for (size_t j = 0; j < xsize - 1; ++j)
          {
            for (size_t k = 0; k < zsize - 1; ++k)
              {
                const size_t index = ModelGeometry.IndexToOffset(j, ysize - 1,
                    k);
                Result(index, index) = -2.0;
                Result(index, ModelGeometry.IndexToOffset(j + 1, ysize - 1, k))
                    = 1.0;
                Result(index, ModelGeometry.IndexToOffset(j, ysize - 1, k + 1))
                    = 1.0;
              }
          }

        //now take care of the border in z-direction
        for (size_t j = 0; j < xsize - 1; ++j)
          {
            for (size_t k = 0; k < ysize - 1; ++k)
              {
                const size_t index = ModelGeometry.IndexToOffset(j, k, zsize
                    - 1);
                Result(index, index) = -2.0;
                Result(index, ModelGeometry.IndexToOffset(j + 1, k, zsize - 1))
                    = 1.0;
                Result(index, ModelGeometry.IndexToOffset(j, k + 1, zsize - 1))
                    = 1.0;
              }
          }
        return Result;
      }
  }
