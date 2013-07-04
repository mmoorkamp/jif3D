//============================================================================
// Name        : test_CellBoundaries.cpp
// Author      : Max Moorkamp
// Version     :
// Copyright   : 2008, MM
//============================================================================

/*!
 * \file test_ThreeDModelBase.cpp
 * This file contains all the tests necessary to verify the functionality of ThreeDModelBase.h
 */
#define BOOST_TEST_MODULE CellBoundaries test
#define BOOST_TEST_MAIN ...
#include <boost/test/included/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <cmath>
#include "../Global/FatalException.h"
#include "../Gravity/ThreeDGravityModel.h"
#include "CellBoundaries.h"

//Test the default state of the object
BOOST_AUTO_TEST_CASE(find_fail_test)
  {
    double Coordinate = 100;
    const size_t ncells = 5;
    jif3D::ThreeDModelBase::t3DModelDim CellBoundaries(boost::extents[ncells]);
    jif3D::ThreeDModelBase::t3DModelDim CellSizes(boost::extents[ncells]);
    std::fill_n(CellSizes.origin(), ncells, 1);
    std::partial_sum(CellSizes.origin(),
        CellSizes.origin() + CellSizes.num_elements() - 1, CellBoundaries.origin() + 1);
    CellBoundaries[0] = 0.0;

    BOOST_CHECK_THROW(jif3D::FindNearestCellBoundary(-1, CellBoundaries, CellSizes),
        jif3D::FatalException);
    BOOST_CHECK_THROW(
        jif3D::FindNearestCellBoundary(Coordinate, CellBoundaries, CellSizes),
        jif3D::FatalException);
  }

BOOST_AUTO_TEST_CASE(find_index_test)
  {

    const size_t ncells = 5;
    jif3D::ThreeDModelBase::t3DModelDim CellBoundaries(boost::extents[ncells]);
    jif3D::ThreeDModelBase::t3DModelDim CellSizes(boost::extents[ncells]);
    std::fill_n(CellSizes.origin(), ncells, 1.0);
    std::partial_sum(CellSizes.origin(),
        CellSizes.origin() + CellSizes.num_elements() - 1, CellBoundaries.origin() + 1);
    CellBoundaries[0] = 0.0;
    size_t index = jif3D::FindNearestCellBoundary(0.1, CellBoundaries, CellSizes);

    size_t nvalues = 50;
    for (size_t i = 0; i < nvalues; ++i)
      {
        //for the last cell the returned index is always the upper bound
        //so we only generate coordinates here that will be rounded
        //to smaller values below and test for the other case separately
        double coordinate = double(i) / double(nvalues) * 4.4;
        size_t index = jif3D::FindNearestCellBoundary(coordinate, CellBoundaries,
            CellSizes);
        BOOST_CHECK_EQUAL(index, std::round(coordinate));
      }
    index = jif3D::FindNearestCellBoundary(4.7, CellBoundaries, CellSizes);
    BOOST_CHECK_EQUAL(index, 4);
  }

BOOST_AUTO_TEST_CASE(index_test)
  {
    jif3D::ThreeDGravityModel Model;
    const size_t nx = 2, ny = 2, nz = 7;
    Model.SetMeshSize(nx, ny, nz);
    std::fill_n(Model.SetXCellSizes().origin(), nx, 1.0);
    std::fill_n(Model.SetYCellSizes().origin(), ny, 1.0);
    std::fill_n(Model.SetZCellSizes().origin(), nz, 1.0);
    size_t nmeas = 70;
    double delta = 0.1;
    for (size_t i = 0; i < nmeas; ++i)
      {
        Model.AddMeasurementPoint(0.5, 0.5, i * delta);
      }
    std::vector<size_t> MeasDepthIndices;
    std::vector<double> ShiftDepth;
    size_t nlevels = jif3D::ConstructDepthIndices(MeasDepthIndices, ShiftDepth, Model);
    BOOST_CHECK_EQUAL(nlevels, nz);
    BOOST_CHECK_EQUAL(MeasDepthIndices.size(), nmeas);
    BOOST_CHECK_EQUAL(ShiftDepth.size(), nz);
    for (size_t i = 0; i < nz; ++i)
      {
        BOOST_CHECK_EQUAL(ShiftDepth[i], i);
      }
    //as above we can calculate the expected values for
    //all measurements above the center of the last cell
    //below that they are always shifted up, so test separately
    for (size_t i = 0; i < nmeas - 5; ++i)
      {
        BOOST_CHECK_EQUAL(MeasDepthIndices[i], std::round(Model.GetMeasPosZ()[i]));
      }
    for (size_t i = nmeas - 5; i < nmeas; ++i)
      {
        BOOST_CHECK_EQUAL(MeasDepthIndices[i], std::floor(Model.GetMeasPosZ()[i]));
      }
  }
