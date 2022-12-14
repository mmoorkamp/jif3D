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
#include "../Gravity/test_common.h"
#include "../Gravity/ThreeDGravityModel.h"
#include "../Global/Jif3DPlatformHelper.h"
#include "EqualGeometry.h"

//Test the default state of the object
BOOST_AUTO_TEST_CASE(constructors_test)
  {
    const jif3D::ThreeDModelBase ConstBaseTest; // 1 //
    BOOST_CHECK_EQUAL(ConstBaseTest.GetXCellSizes().size(), (size_t )0);
    BOOST_CHECK_EQUAL(ConstBaseTest.GetYCellSizes().size(), (size_t )0);
    BOOST_CHECK_EQUAL(ConstBaseTest.GetZCellSizes().size(), (size_t )0);
  }

BOOST_AUTO_TEST_CASE(measpos_test)
  {
    jif3D::ThreeDModelBase BaseTest;
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
        BOOST_CHECK_CLOSE(MeasX.at(i) - oldshiftx, BaseTest.GetMeasPosX().at(i),
            std::numeric_limits<float>::epsilon());
        BOOST_CHECK_CLOSE(MeasY.at(i) - oldshifty, BaseTest.GetMeasPosY().at(i),
            std::numeric_limits<float>::epsilon());
        BOOST_CHECK_CLOSE(MeasZ.at(i) - oldshiftz, BaseTest.GetMeasPosZ().at(i),
            std::numeric_limits<float>::epsilon());
      }
    const double newshiftx = rand();
    const double newshifty = rand();
    const double newshiftz = rand();
    BaseTest.SetOrigin(newshiftx, newshifty, newshiftz);
    for (size_t i = 0; i < nmeas; ++i)
      {
        BOOST_CHECK_CLOSE(MeasX.at(i) - newshiftx, BaseTest.GetMeasPosX().at(i),
            std::numeric_limits<float>::epsilon());
        BOOST_CHECK_CLOSE(MeasY.at(i) - newshifty, BaseTest.GetMeasPosY().at(i),
            std::numeric_limits<float>::epsilon());
        BOOST_CHECK_CLOSE(MeasZ.at(i) - newshiftz, BaseTest.GetMeasPosZ().at(i),
            std::numeric_limits<float>::epsilon());
      }
  }

BOOST_AUTO_TEST_CASE(equal_geometry_test)
  {
    jif3D::ThreeDGravityModel Model1;
    srand((unsigned int)time(nullptr));
    MakeRandomModel(Model1, rand() % 20, 1);
    //we cannot use a ThreeDBaseModel as a concrete object
    //so we use a derived ThreeDGravityModel instead
    jif3D::ThreeDGravityModel Model2(Model1);
    //check that we recognize to equal model geometries as equal
    BOOST_CHECK(jif3D::EqualGridGeometry(Model1, Model2));
    //do we recognize differences in x-direction
    const size_t xsize = Model1.GetXCellSizes().size();
    jif3D::ThreeDGravityModel DiffX(Model1);
    DiffX.SetXCellSizes()[rand() % xsize] += 0.1;
    BOOST_CHECK(!jif3D::EqualGridGeometry(Model1, DiffX));

    //do we recognize differences in y-direction
    const size_t ysize = Model1.GetYCellSizes().size();
    jif3D::ThreeDGravityModel DiffY(Model1);
    DiffY.SetYCellSizes()[rand() % ysize] += 0.1;
    BOOST_CHECK(!jif3D::EqualGridGeometry(Model1, DiffY));

    //do we recognize differences in z-direction
    const size_t zsize = Model1.GetZCellSizes().size();
    jif3D::ThreeDGravityModel DiffZ(Model1);
    DiffZ.SetZCellSizes()[rand() % zsize] += 0.1;
    BOOST_CHECK(!jif3D::EqualGridGeometry(Model1, DiffZ));
  }

BOOST_AUTO_TEST_CASE(concurrent_coordinates_test)
  {
    jif3D::ThreeDGravityModel Model;
    srand((unsigned int) time(nullptr));
    MakeRandomModel(Model, rand() % 20, 1);
    jif3D::ThreeDGravityModel ModelCopy(Model);

#pragma omp parallel default(shared)
      {
#pragma omp for
        for (int i = 0; i < 4; ++i)
          {
            Model.GetXCoordinates();
            Model.GetYCoordinates();
            Model.GetZCoordinates();
          }
      }
    BOOST_CHECK(
        std::equal(Model.GetXCoordinates().begin(), Model.GetXCoordinates().end(),
            ModelCopy.GetXCoordinates().begin()));
    BOOST_CHECK(
        std::equal(Model.GetYCoordinates().begin(), Model.GetYCoordinates().end(),
            ModelCopy.GetYCoordinates().begin()));
    BOOST_CHECK(
        std::equal(Model.GetZCoordinates().begin(), Model.GetZCoordinates().end(),
            ModelCopy.GetZCoordinates().begin()));

  }
