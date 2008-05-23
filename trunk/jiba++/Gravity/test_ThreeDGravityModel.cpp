//============================================================================
// Name        : test_ThreeDGravityModel.cpp
// Author      : Max Moorkamp
// Version     :
// Copyright   : 2008, MM
//============================================================================

/*!
 * \file This file contains all the tests necessary to verify the functionality of ThreeDGravityModel.h
 */
#define BOOST_TEST_MODULE ThreeDGravityModel test
#define BOOST_TEST_MAIN ...
#include <boost/test/included/unit_test.hpp>
#include <stdlib.h>
#include <time.h>
#include <ThreeDGravityModel.h>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/assign/std/vector.hpp> 

using namespace boost::assign;

BOOST_AUTO_TEST_SUITE( Gravity_Test_Suite )
//a helper function
jiba::ThreeDModelBase::t3DModelDim GenerateDimension()
  {

    const size_t DimLength = rand() % 20 +10;
    jiba::ThreeDModelBase::t3DModelDim TestDim(boost::extents[DimLength]);
    for (size_t i = 0; i < DimLength; ++i)
      {
        TestDim[i] = rand() % 25 + 10;
      }
    return TestDim;
  }

//Test the default state of the object
BOOST_AUTO_TEST_CASE(constructors_test)
  {
    const jiba::ThreeDGravityModel ConstBaseTest; // 1 //
    BOOST_CHECK_EQUAL(ConstBaseTest.GetXCellSizes().size(), (size_t)0 );
    BOOST_CHECK_EQUAL(ConstBaseTest.GetYCellSizes().size(), (size_t)0 );
    BOOST_CHECK_EQUAL(ConstBaseTest.GetZCellSizes().size(), (size_t)0 );
    BOOST_CHECK_EQUAL(ConstBaseTest.GetDensities().size(), (size_t)0 );
  }

BOOST_AUTO_TEST_CASE(netcdf_read_write_test)
  {
    srand(time(NULL));
    jiba::ThreeDGravityModel GravityTest;
    jiba::ThreeDModelBase::t3DModelDim XDim = GenerateDimension();
    jiba::ThreeDModelBase::t3DModelDim YDim = GenerateDimension();
    jiba::ThreeDModelBase::t3DModelDim ZDim = GenerateDimension();
    const size_t xsize = XDim.size();
    const size_t ysize = YDim.size();
    const size_t zsize = ZDim.size();
    GravityTest.SetXCellSizes().resize(boost::extents[xsize]);
    GravityTest.SetXCellSizes() = XDim;
    GravityTest.SetYCellSizes().resize(boost::extents[ysize]);
    GravityTest.SetYCellSizes() = YDim;
    GravityTest.SetZCellSizes().resize(boost::extents[zsize]);
    GravityTest.SetZCellSizes() = ZDim;
    BOOST_CHECK(std::equal(XDim.begin(), XDim.end(), GravityTest.GetXCellSizes().begin()));
    BOOST_CHECK(std::equal(YDim.begin(), YDim.end(), GravityTest.GetYCellSizes().begin()));
    BOOST_CHECK(std::equal(ZDim.begin(), ZDim.end(), GravityTest.GetZCellSizes().begin()));

    GravityTest.SetDensities().resize(boost::extents[xsize][ysize][zsize]);
    jiba::ThreeDModelBase::t3DModelData TestData(
        boost::extents[xsize][ysize][zsize]);
    for (size_t i = 0; i < xsize; ++i)
      for (size_t j = 0; j < ysize; ++j)
        for (size_t k = 0; k < zsize; ++k)
          TestData[i][j][k] = rand() % 50;
    GravityTest.SetDensities() = TestData;
    GravityTest.WriteNetCDF("test.nc");
    jiba::ThreeDGravityModel NetCDFReadTest;
    NetCDFReadTest.ReadNetCDF("test.nc");
    BOOST_CHECK(std::equal(XDim.begin(), XDim.end(),
        NetCDFReadTest.GetXCellSizes().begin()));
    BOOST_CHECK(std::equal(YDim.begin(), YDim.end(),
        NetCDFReadTest.GetYCellSizes().begin()));
    BOOST_CHECK(std::equal(ZDim.begin(), ZDim.end(),
        NetCDFReadTest.GetZCellSizes().begin()));
    BOOST_CHECK(std::equal(TestData.begin(), TestData.end(),
        NetCDFReadTest.GetDensities().begin()));
    TestData[0][0][0] += 1.0;
    BOOST_CHECK(!std::equal(TestData.begin(), TestData.end(),
        NetCDFReadTest.GetDensities().begin()));
  }

BOOST_AUTO_TEST_CASE(box_gravity_calc_test)
  {
    //the gravity in the center of the cube should be 0
    BOOST_CHECK(fabs(jiba::CalcGravBox(0.0, 0.0, 0.0, -10.0, -10.0, -10.0,
        20.0, 20.0, 20.0, 1.0)) < std::numeric_limits<double>::epsilon());
    //compare with the reported value of li and chouteau within a precision of 0.1%
    double topofcube = jiba::CalcGravBox(0.0, 0.0, 10.0, -10.0, -10.0, -10.0,
        20.0, 20.0, 20.0, 1.0);
    BOOST_CHECK_CLOSE(topofcube, -3.46426e-6, 0.1);
    double bottomofcube = jiba::CalcGravBox(0.0, 0.0, -10.0, -10.0, -10.0,
        -10.0, 20.0, 20.0, 20.0, 1.0);
    //check symmetry of results
    BOOST_CHECK_CLOSE(topofcube, -bottomofcube,
        std::numeric_limits<float>::epsilon());
    double away1 = jiba::CalcGravBox(1e5, 10.0, 10.0, -10.0, -10.0, -10.0,
        20.0, 20.0, 20.0, 1.0);
    BOOST_CHECK(away1 < std::numeric_limits<double>::epsilon());
    double away2 = jiba::CalcGravBox(100.0, 100.0, 10.0, -10.0, -10.0, -10.0,
        20.0, 20.0, 20.0, 1.0);
    BOOST_CHECK_CLOSE(away2, -1.873178e-9, 0.1);
    //check also for a very large cube
    BOOST_CHECK(fabs(jiba::CalcGravBox(0.0, 0.0, 0.0, -1e6, -1e6, -1e6, 2e6,
        2e6, 2e6, 1.0)) < std::numeric_limits<double>::epsilon());
  }

BOOST_AUTO_TEST_CASE(model_gravity_boxcomp_test)
  {
    double boxtopofcube = jiba::CalcGravBox(11.0, 11.0, 0.0, 0.0, 0.0, 0.0,
        20.0, 20.0, 20.0, 1.0);
    jiba::ThreeDGravityModel GravityTest(true); // store sensitivities
    //create a model of 10x10x10 cells with 2m length in each dimension
    const size_t ncells = 10;
    const double cellsize = 2.0;
    GravityTest.SetXCellSizes().resize(boost::extents[ncells]);
    GravityTest.SetYCellSizes().resize(boost::extents[ncells]);
    GravityTest.SetZCellSizes().resize(boost::extents[ncells]);
    for (size_t i = 0; i < ncells; ++i)
      {
        GravityTest.SetXCellSizes()[i] = cellsize;
        GravityTest.SetYCellSizes()[i] = cellsize;
        GravityTest.SetZCellSizes()[i] = cellsize;
      }
    GravityTest.SetDensities().resize(boost::extents[ncells][ncells][ncells]);
    jiba::rvec DensityVector(ncells*ncells*ncells);
    for (size_t i = 0; i < ncells; ++i)
      for (size_t j = 0; j < ncells; ++j)
        for (size_t k = 0; k < ncells; ++k)
          {
            GravityTest.SetDensities()[i][j][k] = 1.0;
            DensityVector(i*(ncells*ncells)+j*ncells+k) = 1.0;
          }
    GravityTest.AddMeasurementPoint(11.0, 11.0, 0.0);
    jiba::ThreeDGravityModel::tScalarMeasVec
        gravmeas(GravityTest.CalcGravity());
    double gridcube = gravmeas[0];
    GravityTest.WriteNetCDF("cube.nc");
    BOOST_CHECK_CLOSE(boxtopofcube, gridcube,
        std::numeric_limits<float>::epsilon());
    //Check scalar sensitivity calculation
    jiba::rmat ScalarSensitivities(GravityTest.GetScalarSensitivities());
    jiba::rvec SensValues(1);
    SensValues = boost::numeric::ublas::prec_prod(ScalarSensitivities,
        DensityVector);
    BOOST_CHECK_CLOSE(gridcube, SensValues(0),
        std::numeric_limits<float>::epsilon());
  }

BOOST_AUTO_TEST_CASE(background_test)
  {
    jiba::ThreeDGravityModel GravityTest;
    double analytic = 2*M_PI*jiba::Grav_const* (500.0 + 5.0 * 4500.0);
    std::vector<double> bg_dens, bg_thick;
    const size_t nmeas = 10;
    for (size_t i = 0; i < nmeas; ++i)
      GravityTest.AddMeasurementPoint(rand() % 50000 + 2e4, rand() % 50000
          + 2e4, 0.0);
    bg_dens += 1.0,1.0,5.0,5.0;
    bg_thick += 200.0,300.0,3500.0,1000.0;
    GravityTest.SetBackgroundDensities(bg_dens);
    GravityTest.SetBackgroundThicknesses(bg_thick);
    jiba::ThreeDGravityModel::tScalarMeasVec
        gravmeas(GravityTest.CalcGravity());
    for (size_t i = 0; i < nmeas; ++i)
      {
        BOOST_CHECK_CLOSE(analytic, gravmeas[i], 0.01);
        BOOST_CHECK_CLOSE(gravmeas[0], gravmeas[i], 0.01);
      }
  }

BOOST_AUTO_TEST_CASE(model_gravity_1danaly_test)
  {
    jiba::ThreeDGravityModel GravityTest(true);
    //create a model with a random number of cells in each direction
    const size_t nhorcells = rand() % 50 +10;
    const size_t nzcells = 10;
    GravityTest.SetXCellSizes().resize(boost::extents[nhorcells]);
    GravityTest.SetYCellSizes().resize(boost::extents[nhorcells]);
    GravityTest.SetZCellSizes().resize(boost::extents[nzcells]);
    jiba::rvec DensityVector(nhorcells*nhorcells*nzcells + 4); // 4 background layers
    for (size_t i = 0; i < nhorcells; ++i) // set the values of the inner cells
      {
        GravityTest.SetXCellSizes()[i] = rand() % 10000 + 1000;
        GravityTest.SetYCellSizes()[i] = rand() % 10000 + 1000;
      }
    for (size_t i = 0; i < nzcells; ++i)
      {
        GravityTest.SetZCellSizes()[i] = 500;
      }

    GravityTest.SetDensities().resize(boost::extents[nhorcells][nhorcells][nzcells]);
    for (size_t i = 0; i < nhorcells; ++i)
      for (size_t j = 0; j < nhorcells; ++j)
        {
          GravityTest.SetDensities()[i][j][0] = 1.0;
          DensityVector(i*nhorcells*nzcells+j*nzcells) = 1.0;
        }
    for (size_t i = 0; i < nhorcells; ++i)
      for (size_t j = 0; j < nhorcells; ++j)
        for (size_t k = 1; k < nzcells; ++k)
          {
            GravityTest.SetDensities()[i][j][k] = 5.0;
            DensityVector(i*nhorcells*nzcells+j*nzcells+k) = 5.0;
          }

    const size_t nmeas = 10;
    for (size_t i = 0; i < nmeas; ++i)
      GravityTest.AddMeasurementPoint(rand() % 50000 + 2e4, rand() % 50000
          + 2e4, 0.0);
    std::vector<double> bg_dens, bg_thick;
    bg_dens += 1.0, 1.0, 5.0, 5.0;
    bg_thick += 200.0, 300.0, 3500.0, 2000.0;
    
    DensityVector(nhorcells*nhorcells*nzcells) = 1.0;
    DensityVector(nhorcells*nhorcells*nzcells+1) = 1.0;
    DensityVector(nhorcells*nhorcells*nzcells+2) = 5.0;
    DensityVector(nhorcells*nhorcells*nzcells+3) = 5.0;
    GravityTest.SetBackgroundDensities(bg_dens);
    GravityTest.SetBackgroundThicknesses(bg_thick);
    jiba::ThreeDGravityModel::tScalarMeasVec
        scalarmeas(GravityTest.CalcGravity());
    jiba::ThreeDGravityModel::tTensorMeasVec
        tensormeas(GravityTest.CalcTensorGravity());
    GravityTest.WriteNetCDF("layer.nc");
    double analytic = 2*M_PI*jiba::Grav_const* (500.0 + 5.0 * 5500.0);
    jiba::rmat ScalarSensitivities(GravityTest.GetScalarSensitivities());
    jiba::rvec SensValues(prec_prod(ScalarSensitivities, DensityVector));
    //The second time the measurements will also be calculated from the sensitivities internally
    jiba::ThreeDGravityModel::tScalarMeasVec
        scalarmeas2(GravityTest.CalcGravity());
    for (size_t i = 0; i < nmeas; ++i)
      {
        BOOST_CHECK_CLOSE(analytic, scalarmeas[i], 0.01);
        BOOST_CHECK_CLOSE(scalarmeas[i], scalarmeas2[i], 0.001);
        BOOST_CHECK_CLOSE(analytic, SensValues[i], 0.01);
        BOOST_CHECK_CLOSE(scalarmeas[0], scalarmeas[i], 0.01);
        // The  tensorial elements should all be zero
        BOOST_CHECK(fabs(tensormeas[i](0, 0))
            < std::numeric_limits<double>::epsilon());
        BOOST_CHECK(fabs(tensormeas[i](1, 1))
            < std::numeric_limits<double>::epsilon());
        BOOST_CHECK(fabs(tensormeas[i](2, 2))
            < std::numeric_limits<double>::epsilon());
        BOOST_CHECK(fabs(tensormeas[i](0, 1))
            < std::numeric_limits<double>::epsilon());
        BOOST_CHECK(fabs(tensormeas[i](0, 2))
            < std::numeric_limits<double>::epsilon());
        BOOST_CHECK(fabs(tensormeas[i](1, 0))
            < std::numeric_limits<double>::epsilon());
        BOOST_CHECK(fabs(tensormeas[i](1, 2))
            < std::numeric_limits<double>::epsilon());
        BOOST_CHECK(fabs(tensormeas[i](2, 1))
            < std::numeric_limits<double>::epsilon());
        BOOST_CHECK(fabs(tensormeas[i](2, 0))
            < std::numeric_limits<double>::epsilon());
      }
  }

BOOST_AUTO_TEST_CASE(caching_test)
  {
    srand(time(NULL));
    jiba::ThreeDGravityModel GravityTest(true, true);
    jiba::ThreeDModelBase::t3DModelDim XDim = GenerateDimension();
    jiba::ThreeDModelBase::t3DModelDim YDim = GenerateDimension();
    jiba::ThreeDModelBase::t3DModelDim ZDim = GenerateDimension();
    const size_t xsize = XDim.size();
    const size_t ysize = YDim.size();
    const size_t zsize = ZDim.size();
    GravityTest.SetXCellSizes().resize(boost::extents[xsize]);
    GravityTest.SetXCellSizes() = XDim;
    GravityTest.SetYCellSizes().resize(boost::extents[ysize]);
    GravityTest.SetYCellSizes() = YDim;
    GravityTest.SetZCellSizes().resize(boost::extents[zsize]);
    GravityTest.SetZCellSizes() = ZDim;
    std::vector<double> bg_dens, bg_thick;
    bg_dens += 1.0, 1.0, 5.0, 5.0;
    bg_thick += 200.0, 300.0, 3500.0, 2000.0;
    GravityTest.SetBackgroundDensities(bg_dens);
    GravityTest.SetBackgroundThicknesses(bg_thick);
    
    GravityTest.SetDensities().resize(boost::extents[xsize][ysize][zsize]);
    jiba::ThreeDModelBase::t3DModelData TestData(
        boost::extents[xsize][ysize][zsize]);
    for (size_t i = 0; i < xsize; ++i)
      for (size_t j = 0; j < ysize; ++j)
        for (size_t k = 0; k < zsize; ++k)
          TestData[i][j][k] = rand() % 50 + 1;
    GravityTest.SetDensities() = TestData;
    const size_t nmeas = 10;
    for (size_t i = 0; i < nmeas; ++i)
      GravityTest.AddMeasurementPoint(rand() % 100 + 100, rand() % 100 + 100, (rand() % 100 + 100));
    //Calculate twice, once with normal calculation, once cached
    jiba::ThreeDGravityModel::tScalarMeasVec
        scalarmeas1(GravityTest.CalcGravity());
    jiba::ThreeDGravityModel::tScalarMeasVec
        scalarmeas2(GravityTest.CalcGravity());
    for (size_t i = 0; i < nmeas; ++i)
      {
        BOOST_CHECK_CLOSE(scalarmeas1[i],scalarmeas2[i],std::numeric_limits<float>::epsilon());
      }
  }
BOOST_AUTO_TEST_SUITE_END()
