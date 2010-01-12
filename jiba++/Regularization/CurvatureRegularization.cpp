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

    void CurvatureRegularization::ConstructOperator(const jiba::ThreeDModelBase &ModelGeometry)
      {
        const size_t xsize = ModelGeometry.GetModelShape()[0];
                const size_t ysize = ModelGeometry.GetModelShape()[1];
                const size_t zsize = ModelGeometry.GetModelShape()[2];
                //the inner part of the matrix where we can do the derivative in all directions
                for (size_t i = 1; i < xsize - 1; ++i)
                  {
                    for (size_t j = 1; j < ysize - 1; ++j)
                      {
                        for (size_t k = 1; k < zsize - 1; ++k)
                          {
                            const size_t index = ModelGeometry.IndexToOffset(i, j, k);
                            XOperatorMatrix(index, ModelGeometry.IndexToOffset(i + 1,
                                j, k)) = 1.0;
                            XOperatorMatrix(index, index) = -2.0;
                            XOperatorMatrix(index, ModelGeometry.IndexToOffset(i - 1,
                                j, k)) = 1.0;

                            YOperatorMatrix(index, ModelGeometry.IndexToOffset(i,
                                j + 1, k)) = 1.0;
                            YOperatorMatrix(index, index) = -2.0;
                            YOperatorMatrix(index, ModelGeometry.IndexToOffset(i,
                                j - 1, k)) = 1.0;

                            ZOperatorMatrix(index, ModelGeometry.IndexToOffset(i, j, k
                                + 1)) = 1.0;
                            ZOperatorMatrix(index, index) = -2.0;
                            ZOperatorMatrix(index, ModelGeometry.IndexToOffset(i, j, k
                                - 1)) = 1.0;
                          }
                        //we can handle the border in z-direction within this loop
                        //first the maximum z index
                        size_t index = ModelGeometry.IndexToOffset(i, j, zsize - 1);
                        XOperatorMatrix(index, ModelGeometry.IndexToOffset(i + 1, j,
                            zsize - 1)) = 1.0;
                        XOperatorMatrix(index, index) = -2.0;
                        XOperatorMatrix(index, ModelGeometry.IndexToOffset(i - 1, j,
                            zsize - 1)) = 1.0;

                        YOperatorMatrix(index, ModelGeometry.IndexToOffset(i, j + 1,
                            zsize - 1)) = 1.0;
                        YOperatorMatrix(index, index) = -2.0;
                        YOperatorMatrix(index, ModelGeometry.IndexToOffset(i, j - 1,
                            zsize - 1)) = 1.0;

                        //then the first z-index
                        index = ModelGeometry.IndexToOffset(i, j, 0);
                        XOperatorMatrix(index, ModelGeometry.IndexToOffset(i + 1, j, 0))
                            = 1.0;
                        XOperatorMatrix(index, index) = -2.0;
                        XOperatorMatrix(index, ModelGeometry.IndexToOffset(i - 1, j, 0))
                            = 1.0;

                        YOperatorMatrix(index, ModelGeometry.IndexToOffset(i, j + 1, 0))
                            = 1.0;
                        YOperatorMatrix(index, index) = -2.0;
                        YOperatorMatrix(index, ModelGeometry.IndexToOffset(i, j - 1, 0))
                            = 1.0;
                      }
                  }
      }
  }
