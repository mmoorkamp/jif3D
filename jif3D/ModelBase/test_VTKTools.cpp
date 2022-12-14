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
#include <boost/test/tools/floating_point_comparison.hpp>
#include <numeric>
#include <cstdlib>
#include "../Gravity/ThreeDGravityModel.h"
#include "../Global/NumUtil.h"
#include "../Global/Jif3DPlatformHelper.h"
#include "VTKTools.h"

//Test the default state of the object
BOOST_AUTO_TEST_CASE(read_write_test)
  {
    srand((unsigned int) time(0));
    jif3D::ThreeDGravityModel Model, Compare;
    size_t nx = rand() % 30;
    size_t ny = rand() % 30;
    size_t nz = rand() % 30;
    jif3D::ThreeDModelBase::t3DModelDim XCS(nx), YCS(ny), ZCS(nz);

    Model.SetDensities().resize(boost::extents[nx][ny][nz]);
    std::generate_n(XCS.begin(), nx, jif3D::platform::drand48);
    std::generate_n(YCS.begin(), ny, jif3D::platform::drand48);
    std::generate_n(ZCS.begin(), nz, jif3D::platform::drand48);
    Model.SetXCellSizes(XCS);
    Model.SetYCellSizes(YCS);
    Model.SetZCellSizes(ZCS);
    std::generate_n(Model.SetDensities().origin(), nx * ny * nz,
        jif3D::platform::drand48);
    std::string filename = "rwtest.vtk";
    jif3D::Write3DModelToVTK(filename, "Density", Model.GetXCoordinates(),
        Model.GetYCoordinates(), Model.GetZCoordinates(), Model.GetDensities());
    jif3D::Read3DModelFromVTK(filename, XCS, YCS, ZCS, Compare.SetDensities());
    Compare.SetXCoordinates(XCS);
    Compare.SetYCoordinates(YCS);
    Compare.SetZCoordinates(ZCS);
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

