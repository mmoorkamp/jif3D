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

    GradientRegularization::comp_mat GradientRegularization::SmoothnessOperator()
      {
        comp_mat Result(ModelGeometry.GetData().num_elements(),
            ModelGeometry.GetData().num_elements());
        const size_t xsize = ModelGeometry.GetData().shape()[0];
        const size_t ysize = ModelGeometry.GetData().shape()[1];
        const size_t zsize = ModelGeometry.GetData().shape()[2];
        for (size_t i = 1; i < xsize - 1; ++i)
          {
            for (size_t j = 1; j < ysize - 1; ++j)
              {
                for (size_t k = 1; k < zsize - 1; ++k)
                  {
                    const size_t index = ModelGeometry.IndexToOffset(i, j, k);
                    Result(index, ModelGeometry.IndexToOffset(i + 1, j, k))
                        = 1.0;
                    Result(index, ModelGeometry.IndexToOffset(i - 1, j, k))
                        = -1.0;

                    Result(index, ModelGeometry.IndexToOffset(i, j + 1, k))
                        = 1.0;
                    Result(index, ModelGeometry.IndexToOffset(i, j - 1, k))
                        = -1.0;

                    Result(index, ModelGeometry.IndexToOffset(i, j, k + 1))
                        = 1.0;
                    Result(index, ModelGeometry.IndexToOffset(i, j, k - 1))
                        = -1.0;

                  }
              }
          }
        return Result;
      }
  }
