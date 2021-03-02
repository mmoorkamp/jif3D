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
#include "../Global/Jif3DTesting.h"
#include <boost/test/included/unit_test.hpp>
#include <boost/test/tools/floating_point_comparison.hpp>
#include <numeric>
#include "../Gravity/test_common.h"
#include "../Gravity/ThreeDGravityModel.h"
#include "../Gravity/ScalarGravityData.h"
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

BOOST_AUTO_TEST_CASE(equal_geometry_test)
  {
    jif3D::ThreeDGravityModel Model1;
    jif3D::ScalarGravityData Data;
    srand((unsigned int) time(nullptr));
    MakeRandomModel(Model1, Data, rand() % 20, 1);
    //we cannot use a ThreeDBaseModel as a concrete object
    //so we use a derived ThreeDGravityModel instead
    jif3D::ThreeDGravityModel Model2(Model1);
    //check that we recognize to equal model geometries as equal
    BOOST_CHECK(jif3D::EqualGridGeometry(Model1, Model2));
    //do we recognize differences in x-direction
    const size_t xsize = Model1.GetXCellSizes().size();
    jif3D::ThreeDGravityModel DiffX(Model1);
    auto XCS = Model1.GetXCellSizes();
    XCS[rand() % xsize] += 0.1;
    DiffX.SetXCellSizes(XCS);
    BOOST_CHECK(!jif3D::EqualGridGeometry(Model1, DiffX));

    //do we recognize differences in y-direction
    const size_t ysize = Model1.GetYCellSizes().size();
    jif3D::ThreeDGravityModel DiffY(Model1);
    auto YCS = Model1.GetYCellSizes();
    YCS[rand() % ysize] += 0.1;
    DiffY.SetYCellSizes(YCS);
    BOOST_CHECK(!jif3D::EqualGridGeometry(Model1, DiffY));

    //do we recognize differences in z-direction
    const size_t zsize = Model1.GetZCellSizes().size();
    jif3D::ThreeDGravityModel DiffZ(Model1);
    auto ZCS = Model1.GetZCellSizes();
    ZCS[rand() % zsize] += 0.1;
    DiffZ.SetZCellSizes(ZCS);
    BOOST_CHECK(!jif3D::EqualGridGeometry(Model1, DiffZ));

    //check differences in x-origin
    jif3D::ThreeDGravityModel OriginX(Model1);
    auto XCoord = Model1.GetXCoordinates();
    std::transform(XCoord.begin(), XCoord.end(), XCoord.begin(), [](double val)
      { return val + 1.0;});
    OriginX.SetXCoordinates(XCoord);
    BOOST_CHECK(!jif3D::EqualGridGeometry(Model1, OriginX));

    //check differences in y-origin
    jif3D::ThreeDGravityModel OriginY(Model1);
    auto YCoord = Model1.GetYCoordinates();
    std::transform(YCoord.begin(), YCoord.end(), YCoord.begin(), [](double val)
      { return val + 2.0;});
    OriginY.SetYCoordinates(YCoord);
    BOOST_CHECK(!jif3D::EqualGridGeometry(Model1, OriginY));

    //check differences in origin
    jif3D::ThreeDGravityModel OriginZ(Model1);
    auto ZCoord = Model1.GetZCoordinates();
    std::transform(ZCoord.begin(), ZCoord.end(), ZCoord.begin(), [](double val)
      { return val + 3.0;});
    OriginZ.SetZCoordinates(ZCoord);
    BOOST_CHECK(!jif3D::EqualGridGeometry(Model1, OriginZ));
  }

BOOST_AUTO_TEST_CASE(concurrent_coordinates_test)
  {
    jif3D::ThreeDGravityModel Model;
    jif3D::ScalarGravityData Data;
    srand((unsigned int) time(nullptr));
    MakeRandomModel(Model, Data, rand() % 20, 1);
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
