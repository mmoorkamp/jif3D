//============================================================================
// Name        : test_ModelRefiner.cpp
// Author      : Feb 20, 2009
// Version     :
// Copyright   : 2009, mmoorkamp
//============================================================================


#define BOOST_TEST_MODULE ModelRefiner test
#define BOOST_TEST_MAIN ...
#include <boost/test/included/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include "ModelRefiner.h"
#include "ThreeDModelBase.h"
#include "../Gravity/ThreeDGravityModel.h"

BOOST_AUTO_TEST_CASE(basic_axis_refinement)
  {
    jiba::ThreeDGravityModel CoarseModel;
    jiba::ThreeDGravityModel FineModel;
    size_t ncoarsecells = 2;
    CoarseModel.SetXCellSizes().resize(boost::extents[ncoarsecells]);
    CoarseModel.SetYCellSizes().resize(boost::extents[ncoarsecells]);
    CoarseModel.SetZCellSizes().resize(boost::extents[ncoarsecells]);
    for (size_t i = 0; i < ncoarsecells; ++i)
      {
        CoarseModel.SetXCellSizes()[i] = 50;
        CoarseModel.SetYCellSizes()[i] = 60;
        CoarseModel.SetZCellSizes()[i] = 70;
      }

    const size_t xrefinement = 10;
    jiba::ThreeDModelBase::t3DModelDim XRefiner(boost::extents[xrefinement]);
    for (size_t i = 0; i < xrefinement; ++i)
      {
        XRefiner[i] = 17 * i;
      }
    jiba::ModelRefiner Refiner;
    Refiner.SetXCoordinates(XRefiner);
    Refiner.RefineAxes(CoarseModel,FineModel);
    for (size_t i = 0; i < FineModel.GetXCellSizes().size(); ++i)
      std::cout << FineModel.GetXCellSizes()[i] << std::endl;
}
