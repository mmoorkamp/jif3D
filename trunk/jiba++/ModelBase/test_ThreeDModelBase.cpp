//============================================================================
// Name        : test_ThreeDModelBase.cpp
// Author      : Max Moorkamp
// Version     :
// Copyright   : 2008, MM
//============================================================================

/*!
 * \file This file contains all the tests necessary to verify the functionality of ThreeDModelBase.h
 */
#define BOOST_TEST_MODULE ThreeDModelBase test
#define BOOST_TEST_MAIN ...
#include <boost/test/included/unit_test.hpp>

#include <ThreeDModelBase.h>

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
