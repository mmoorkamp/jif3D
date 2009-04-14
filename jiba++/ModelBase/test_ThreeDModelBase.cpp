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
#include <boost/test/floating_point_comparison.hpp>
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


BOOST_AUTO_TEST_CASE(measpos_test)
  {
    jiba::ThreeDModelBase BaseTest;
    const double oldshiftx = rand();
    const double oldshifty = rand();
    const double oldshiftz = rand();
    const size_t nmeas = 10;
    std::vector<double> MeasX, MeasY, MeasZ;
    for (size_t i = 0; i < nmeas; ++i)
      {
        MeasX.push_back(rand());
        MeasY.push_back(rand());
        MeasZ.push_back(rand());
        BaseTest.AddMeasurementPoint(MeasX.at(i), MeasY.at(i), MeasZ.at(i));
      }
    BaseTest.SetOrigin(oldshiftx, oldshifty, oldshiftz);
    for (size_t i = 0; i < nmeas; ++i)
      {
        BOOST_CHECK_CLOSE (MeasX.at(i) - oldshiftx, BaseTest.GetMeasPosX().at(i), std::numeric_limits<
            float>::epsilon());
        BOOST_CHECK_CLOSE (MeasY.at(i) - oldshifty, BaseTest.GetMeasPosY().at(i), std::numeric_limits<
            float>::epsilon());
        BOOST_CHECK_CLOSE (MeasZ.at(i) - oldshiftz, BaseTest.GetMeasPosZ().at(i), std::numeric_limits<
            float>::epsilon());
      }
    const double newshiftx = rand();
    const double newshifty = rand();
    const double newshiftz = rand();
    BaseTest.SetOrigin(newshiftx, newshifty, newshiftz);
    for (size_t i = 0; i < nmeas; ++i)
      {
        BOOST_CHECK_CLOSE (MeasX.at(i) - newshiftx, BaseTest.GetMeasPosX().at(i), std::numeric_limits<
            float>::epsilon());
        BOOST_CHECK_CLOSE (MeasY.at(i) - newshifty, BaseTest.GetMeasPosY().at(i), std::numeric_limits<
            float>::epsilon());
        BOOST_CHECK_CLOSE (MeasZ.at(i) - newshiftz, BaseTest.GetMeasPosZ().at(i), std::numeric_limits<
            float>::epsilon());
      }
  }
