//============================================================================
// Name        : test_ThreeDModelBase.cpp
// Author      : Max Moorkamp
// Version     :
// Copyright   : 2008, MM
//============================================================================

/*!
 * \file test_ThreeDModelBase.cpp
 * This file contains all the tests necessary to verify the functionality of ThreeDModelBase.h
 */
#define BOOST_TEST_MODULE ThreeDModelBase test
#define BOOST_TEST_MAIN ...
#include <boost/test/included/unit_test.hpp>
#include <numeric>
#include "ThreeDModelBase.h"

//Test the default state of the object
BOOST_AUTO_TEST_CASE(constructors_test)
  {
    const jiba::ThreeDModelBase ConstBaseTest; // 1 //
    BOOST_CHECK_EQUAL(ConstBaseTest.GetXCellSizes().size(), (size_t)0 );
    BOOST_CHECK_EQUAL(ConstBaseTest.GetYCellSizes().size(), (size_t)0 );
    BOOST_CHECK_EQUAL(ConstBaseTest.GetZCellSizes().size(), (size_t)0 );
  }

//Test the access function for the cell dimensions
BOOST_AUTO_TEST_CASE(cell_dimension_access_function_test)
  {
    jiba::ThreeDModelBase BaseTest;
    jiba::ThreeDModelBase::t3DModelDim TestDim(boost::extents[5]);
    for (size_t i = 0; i < 5; ++i)
      TestDim[i] = i;
    BaseTest.SetXCellSizes().resize(boost::extents[5]);
    BaseTest.SetXCellSizes() = TestDim;
    BaseTest.SetYCellSizes().resize(boost::extents[5]);
    BaseTest.SetYCellSizes() = TestDim;
    BaseTest.SetZCellSizes().resize(boost::extents[5]);
    BaseTest.SetZCellSizes() = TestDim;
    BOOST_CHECK(std::equal(TestDim.begin(), TestDim.end(),
            BaseTest.GetXCellSizes().begin()));
    BOOST_CHECK(std::equal(TestDim.begin(), TestDim.end(),
            BaseTest.GetYCellSizes().begin()));
    BOOST_CHECK(std::equal(TestDim.begin(), TestDim.end(),
            BaseTest.GetZCellSizes().begin()));
  }

BOOST_AUTO_TEST_CASE(coordinate_calculation_test)
  {
    const size_t nelements = 5;
    jiba::ThreeDModelBase BaseTest;
    jiba::ThreeDModelBase::t3DModelDim TestDim(boost::extents[nelements]);
    for (size_t i = 0; i < nelements; ++i)
      TestDim[i] = rand();
    BaseTest.SetXCellSizes().resize(boost::extents[nelements]);
    BaseTest.SetXCellSizes() = TestDim;
    BaseTest.SetYCellSizes().resize(boost::extents[nelements]);
    BaseTest.SetYCellSizes() = TestDim;
    BaseTest.SetZCellSizes().resize(boost::extents[nelements]);
    BaseTest.SetZCellSizes() = TestDim;
    jiba::ThreeDModelBase::t3DModelDim XCoordinates(BaseTest.GetXCoordinates());
    jiba::ThreeDModelBase::t3DModelDim YCoordinates(BaseTest.GetYCoordinates());
    jiba::ThreeDModelBase::t3DModelDim ZCoordinates(BaseTest.GetZCoordinates());
    BOOST_CHECK_EQUAL(XCoordinates[0],0.0);
    BOOST_CHECK_EQUAL(YCoordinates[0],0.0);
    BOOST_CHECK_EQUAL(ZCoordinates[0],0.0);
    jiba::ThreeDModelBase::t3DModelDim CompCoord(boost::extents[nelements]);
    std::partial_sum(TestDim.begin(), TestDim.end(), CompCoord.begin());
    BOOST_CHECK(std::equal(CompCoord.begin(), CompCoord.end()-1,
            BaseTest.GetXCoordinates().begin()+1));
    BOOST_CHECK(std::equal(CompCoord.begin(), CompCoord.end()-1,
            BaseTest.GetYCoordinates().begin()+1));
    BOOST_CHECK(std::equal(CompCoord.begin(), CompCoord.end()-1,
            BaseTest.GetZCoordinates().begin()+1));
  }

BOOST_AUTO_TEST_CASE(find_indices_test)
  {
    const size_t nelements = 5;
    jiba::ThreeDModelBase BaseTest;
    BaseTest.SetXCellSizes().resize(boost::extents[nelements]);
    BaseTest.SetYCellSizes().resize(boost::extents[nelements]);
    BaseTest.SetZCellSizes().resize(boost::extents[nelements]);
    for (size_t i = 0; i < nelements; ++i)
      {
        BaseTest.SetXCellSizes()[i] = 1.5;
        BaseTest.SetYCellSizes()[i] = 2.0;
        BaseTest.SetZCellSizes()[i] = 2.5;

      }
    boost::array<jiba::ThreeDModelBase::t3DModelData::index,3> indices(BaseTest.FindAssociatedIndices(2.2,2.2,2.2));
    BOOST_CHECK_EQUAL(indices[0], 1);
    BOOST_CHECK_EQUAL(indices[1], 1);
    BOOST_CHECK_EQUAL(indices[2], 0);
  }
