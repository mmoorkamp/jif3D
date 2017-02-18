//============================================================================
// Name        : CurvatureRegularization.cpp
// Author      : Jan 11, 2010
// Version     : 
// Copyright   : 2010, mmoorkamp
//============================================================================

#include "CurvatureRegularization.h"

namespace jif3D
  {

    CurvatureRegularization::~CurvatureRegularization()
      {

      }

    void CurvatureRegularization::ConstructOperator(
        const jif3D::ThreeDModelBase &ModelGeometry,
        const jif3D::ThreeDModelBase &TearModelX,
        const jif3D::ThreeDModelBase &TearModelY,
        const jif3D::ThreeDModelBase &TearModelZ)
      {
        const size_t xsize = ModelGeometry.GetModelShape()[0];
        const size_t ysize = ModelGeometry.GetModelShape()[1];
        const size_t zsize = ModelGeometry.GetModelShape()[2];
        const double CenterValue = -2.0 + Eps;
        const double LeftValue = 1.0;
        const double RightValue = 1.0;
        const double AvgX = std::accumulate(ModelGeometry.GetXCellSizes().begin(),
            ModelGeometry.GetXCellSizes().end(), 0.0)
            / boost::numeric_cast<double>(ModelGeometry.GetXCellSizes().size());
        const double AvgY = std::accumulate(ModelGeometry.GetYCellSizes().begin(),
            ModelGeometry.GetYCellSizes().end(), 0.0)
            / boost::numeric_cast<double>(ModelGeometry.GetYCellSizes().size());
        const double AvgZ = std::accumulate(ModelGeometry.GetZCellSizes().begin(),
            ModelGeometry.GetZCellSizes().end(), 0.0)
            / boost::numeric_cast<double>(ModelGeometry.GetZCellSizes().size());
        //the inner part of the matrix where we can do the derivative in all directions
        for (size_t i = 1; i < xsize - 1; ++i)
          {
            const double XCSize =
                WithCellsize ?
                    std::pow(
                        (ModelGeometry.GetXCellSizes()[i + 1]
                            + ModelGeometry.GetXCellSizes()[i]) / (2.0 * AvgX), 2) :
                    1.0;
            for (size_t j = 0; j < ysize; ++j)
              {
                for (size_t k = 0; k < zsize; ++k)
                  {
                    //save the storage index of the current element in the inversion domain
                    const size_t index = ModelGeometry.IndexToOffset(i, j, k);
                    //for each direction we have to check whether we want to introduce a tear
                    //we have to check the current element and the previous one,
                    //as curvature based regularization depends on elements on both
                    //sides of the current element, if both are different from zero
                    //we add the appropriate values to the operator matrix
                    if (TearModelX.GetData()[i][j][k]
                        && TearModelX.GetData()[i - 1][j][k])
                      {
                        XOperatorMatrix(index, ModelGeometry.IndexToOffset(i + 1, j, k)) =
                            RightValue * std::abs(TearModelX.GetData()[i][j][k]) / XCSize;
                        XOperatorMatrix(index, index) = CenterValue
                            * std::abs(TearModelX.GetData()[i][j][k]) / XCSize;
                        XOperatorMatrix(index, ModelGeometry.IndexToOffset(i - 1, j, k)) =
                            LeftValue * std::abs(TearModelX.GetData()[i][j][k]) / XCSize;
                      }
                  }
              } //end of loop in y-direction
          } //end of loop in z-direction

        for (size_t i = 0; i < xsize; ++i)
          {
            for (size_t j = 1; j < ysize - 1; ++j)
              {
                const double YCSize =
                    WithCellsize ?
                        std::pow(
                            (ModelGeometry.GetYCellSizes()[j + 1]
                                + ModelGeometry.GetYCellSizes()[j]) / (2.0 * AvgY), 2) :
                        1.0;
                for (size_t k = 0; k < zsize; ++k)
                  {
                    //save the storage index of the current element in the inversion domain
                    const size_t index = ModelGeometry.IndexToOffset(i, j, k);
                    //for each direction we have to check whether we want to introduce a tear
                    //we have to check the current element and the previous one,
                    //as curvature based regularization depends on elements on both
                    //sides of the current element, if both are different from zero
                    //we add the appropriate values to the operator matrix
                    if (TearModelY.GetData()[i][j][k]
                        && TearModelY.GetData()[i][j - 1][k])
                      {
                        YOperatorMatrix(index, ModelGeometry.IndexToOffset(i, j + 1, k)) =
                            RightValue * std::abs(TearModelY.GetData()[i][j][k]) / YCSize;
                        YOperatorMatrix(index, index) = CenterValue
                            * std::abs(TearModelY.GetData()[i][j][k]) / YCSize;
                        YOperatorMatrix(index, ModelGeometry.IndexToOffset(i, j - 1, k)) =
                            LeftValue * std::abs(TearModelY.GetData()[i][j][k]) / YCSize;
                      }
                  }
              } //end of loop in y-direction
          } //end of loop in z-direction

        for (size_t i = 0; i < xsize; ++i)
          {
            for (size_t j = 0; j < ysize; ++j)
              {
                for (size_t k = 1; k < zsize - 1; ++k)
                  {
                    const double ZCSize =
                        WithCellsize ?
                            std::pow(
                                (ModelGeometry.GetZCellSizes()[k + 1]
                                    + ModelGeometry.GetZCellSizes()[k]) / (2.0 * AvgZ),
                                2) :
                            1.0;
                    //save the storage index of the current element in the inversion domain
                    const size_t index = ModelGeometry.IndexToOffset(i, j, k);
                    //for each direction we have to check whether we want to introduce a tear
                    //we have to check the current element and the previous one,
                    //as curvature based regularization depends on elements on both
                    //sides of the current element, if both are different from zero
                    //we add the appropriate values to the operator matrix

                    if (TearModelZ.GetData()[i][j][k]
                        && TearModelZ.GetData()[i][j][k - 1])
                      {
                        ZOperatorMatrix(index, ModelGeometry.IndexToOffset(i, j, k + 1)) =
                            RightValue * std::abs(TearModelZ.GetData()[i][j][k]) / ZCSize;
                        ZOperatorMatrix(index, index) = CenterValue
                            * std::abs(TearModelZ.GetData()[i][j][k]) / ZCSize;
                        ZOperatorMatrix(index, ModelGeometry.IndexToOffset(i, j, k - 1)) =
                            LeftValue * std::abs(TearModelZ.GetData()[i][j][k]) / ZCSize;
                      }
                  }
              } //end of loop in y-direction
          } //end of loop in z-direction
      }
  }
