//============================================================================
// Name        : test_VTKTools.cpp
// Author      : Max Moorkamp
// Version     :
// Copyright   : 2008, MM
//============================================================================

/*!
 * \file test_ThreeDModelBase.cpp
 * This file contains all the tests necessary to verify the functionality of ThreeDModelBase.h
 */
#define BOOST_TEST_MODULE VTKTools test
#define BOOST_TEST_MAIN ...
#include "../Global/Jif3DTesting.h"
#include <boost/test/included/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <numeric>
#include <cstdlib>
#include "../Gravity/ThreeDGravityModel.h"
#include "../Global/NumUtil.h"
#include "../Global/Jif3DPlatformHelper.h"
#include "VTKTools.h"

//Test the default state of the object
BOOST_AUTO_TEST_CASE(read_write_test)
  {
    srand((unsigned int)time(0));
    jif3D::ThreeDGravityModel Model, Compare;
    size_t nx = rand() % 30;
    size_t ny = rand() % 30;
    size_t nz = rand() % 30;
    Model.SetXCellSizes().resize(nx);
    Model.SetYCellSizes().resize(ny);
    Model.SetZCellSizes().resize(nz);
    Model.SetDensities().resize(boost::extents[nx][ny][nz]);
    std::generate_n(Model.SetXCellSizes().begin(), nx, jif3D::platform::drand48);
    std::generate_n(Model.SetYCellSizes().begin(), ny, jif3D::platform::drand48);
    std::generate_n(Model.SetZCellSizes().begin(), nz, jif3D::platform::drand48);
    std::generate_n(Model.SetDensities().origin(), nx * ny * nz, jif3D::platform::drand48);
    std::string filename = "rwtest.vtk";
    jif3D::Write3DModelToVTK(filename, "Density", Model.GetXCellSizes(),
        Model.GetYCellSizes(), Model.GetZCellSizes(), Model.GetDensities());
    jif3D::Read3DModelFromVTK(filename, Compare.SetXCellSizes(), Compare.SetYCellSizes(),
        Compare.SetZCellSizes(), Compare.SetDensities());
    BOOST_CHECK_EQUAL(Model.GetXCellSizes().size(), Compare.GetXCellSizes().size());
    BOOST_CHECK_EQUAL(Model.GetYCellSizes().size(), Compare.GetYCellSizes().size());
    BOOST_CHECK_EQUAL(Model.GetZCellSizes().size(), Compare.GetZCellSizes().size());
    BOOST_CHECK_EQUAL(Model.GetDensities().num_elements(),
        Compare.GetDensities().num_elements());
    BOOST_CHECK(
        std::equal(Model.GetXCellSizes().begin(), Model.GetXCellSizes().end(),
            Compare.GetXCellSizes().begin(), jif3D::roughlyEqual<double, double>(1e-3)));
    BOOST_CHECK(
        std::equal(Model.GetYCellSizes().begin(), Model.GetYCellSizes().end(),
            Compare.GetYCellSizes().begin(), jif3D::roughlyEqual<double, double>(1e-3)));
    BOOST_CHECK(
        std::equal(Model.GetZCellSizes().begin(), Model.GetZCellSizes().end(),
            Compare.GetZCellSizes().begin(), jif3D::roughlyEqual<double, double>(1e-3)));
    BOOST_CHECK(
        std::equal(Model.GetDensities().origin(),
            Model.GetDensities().origin() + Model.GetDensities().num_elements(),
            Compare.GetDensities().origin(), jif3D::roughlyEqual<double, double>(1e-3)));
  }

