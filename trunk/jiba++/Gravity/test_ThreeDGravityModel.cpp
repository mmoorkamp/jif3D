//============================================================================
// Name        : test_ThreeDGravityModel.cpp
// Author      : Max Moorkamp
// Version     :
// Copyright   : 2008, MM
//============================================================================

/*!
 * \file test_ThreeDGravityModel.cpp
 * This file contains all the tests necessary to verify the functionality of ThreeDGravityModel.h
 */
#define BOOST_TEST_MODULE ThreeDGravityModel test
#define BOOST_TEST_MAIN ...
#include <boost/test/included/unit_test.hpp>
#include <stdlib.h>
#include <time.h>
#include <fstream>
#include "ThreeDGravityModel.h"
#include "BasicGravElements.h"
#include "MinMemGravityCalculator.h"
#include <boost/test/floating_point_comparison.hpp>
#include <boost/assign/std/vector.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/numeric/conversion/cast.hpp>
#include "ScalarOMPGravityImp.h"

using namespace boost::assign;

BOOST_AUTO_TEST_SUITE( Gravity_Test_Suite )

//a helper function to create a model dimension of random size
jiba  ::ThreeDModelBase::t3DModelDim GenerateDimension()
    {

      const size_t DimLength = rand() % 11 +15;
      jiba::ThreeDModelBase::t3DModelDim TestDim(boost::extents[DimLength]);
      for (size_t i = 0; i < DimLength; ++i)
        {
          TestDim[i] = rand() % 25 + 10;
        }
      return TestDim;
    }

  //create a random density model without any background layers
  //this is necessary to test the R interface, which does not allow
  //for background layers, yet
  void MakeRandomModel(jiba::ThreeDGravityModel &Model, const size_t nmeas)
    {
      srand(time(NULL));

      jiba::ThreeDModelBase::t3DModelDim XDim = GenerateDimension();
      jiba::ThreeDModelBase::t3DModelDim YDim = GenerateDimension();
      jiba::ThreeDModelBase::t3DModelDim ZDim = GenerateDimension();
      const size_t xsize = XDim.size();
      const size_t ysize = YDim.size();
      const size_t zsize = ZDim.size();
      Model.SetXCellSizes().resize(boost::extents[xsize]);
      Model.SetXCellSizes() = XDim;
      Model.SetYCellSizes().resize(boost::extents[ysize]);
      Model.SetYCellSizes() = YDim;
      Model.SetZCellSizes().resize(boost::extents[zsize]);
      Model.SetZCellSizes() = ZDim;
      int xlength = boost::numeric_cast<int>(floor(std::accumulate(XDim.begin(), XDim.end(), 0.0)));
      int ylength = boost::numeric_cast<int>(floor(std::accumulate(YDim.begin(), YDim.end(), 0.0)));

      Model.SetDensities().resize(boost::extents[xsize][ysize][zsize]);

      for (size_t i = 0; i < xsize; ++i)
      for (size_t j = 0; j < ysize; ++j)
      for (size_t k = 0; k < zsize; ++k)
        {
          if (i < xsize / 2)
          Model.SetDensities()[i][j][k] = double(rand() % 50) / 10.0 + 1.0;
          else
          Model.SetDensities()[i][j][k] = -double(rand() % 50) / 10.0 + 1.0;
        }
      //generate measurement  points
      //the z-axis is positive down, so we choose negative z-coordinates
      // => we are measuring above the surface
      for (size_t i = 0; i < nmeas; ++i)
      Model.AddMeasurementPoint(rand() % xlength + 1, rand() % ylength + 1,
          -(rand() % 100 + 1.0));
    }

  //Test the default state of the object
  BOOST_AUTO_TEST_CASE(constructors_test)
    {
      const jiba::ThreeDGravityModel ConstBaseTest; // 1 //
      BOOST_CHECK_EQUAL(ConstBaseTest.GetXCellSizes().size(), (size_t) 0);
      BOOST_CHECK_EQUAL(ConstBaseTest.GetYCellSizes().size(), (size_t) 0);
      BOOST_CHECK_EQUAL(ConstBaseTest.GetZCellSizes().size(), (size_t) 0);
      BOOST_CHECK_EQUAL(ConstBaseTest.GetDensities().size(), (size_t) 0);
    }

  BOOST_AUTO_TEST_CASE(indexoffset_test)
    {
      jiba::ThreeDGravityModel ModelTest;
      size_t xsize = rand() % 100 + 1;
      size_t ysize = rand() % 100 + 1;
      size_t zsize = rand() % 100 +1;
      ModelTest.SetDensities().resize(boost::extents[xsize][ysize][zsize]);
      for (int i = 0; i < 10; ++i)
        {
          size_t xindex = rand() % xsize;
          size_t yindex = rand() % ysize;
          size_t zindex = rand() % zsize;
          int offset = ModelTest.IndexToOffset(xindex,yindex,zindex);
          int newx, newy, newz;
          ModelTest.OffsetToIndex(offset,newx,newy,newz);
          BOOST_CHECK_EQUAL(xindex,newx);
          BOOST_CHECK_EQUAL(yindex,newy);
          BOOST_CHECK_EQUAL(zindex,newz);
        }
    }

  //Check whether storing the model in a netcdf file and reading it in again creates the same model
  BOOST_AUTO_TEST_CASE(netcdf_read_write_test)
    {
      jiba::ThreeDGravityModel GravityTest;

      srand(time(NULL));

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

      GravityTest.SetDensities().resize(boost::extents[xsize][ysize][zsize]);

      jiba::ThreeDModelBase::t3DModelData TestData(
          boost::extents[xsize][ysize][zsize]);
      for (size_t i = 0; i < xsize; ++i)
      for (size_t j = 0; j < ysize; ++j)
      for (size_t k = 0; k < zsize; ++k)
      TestData[i][j][k] = (double(rand() % 50) / 10.0 + 1.0);
      GravityTest.SetDensities() = TestData;
      BOOST_CHECK(std::equal(XDim.begin(), XDim.end(),
              GravityTest.GetXCellSizes().begin()));
      BOOST_CHECK(std::equal(YDim.begin(), YDim.end(),
              GravityTest.GetYCellSizes().begin()));
      BOOST_CHECK(std::equal(ZDim.begin(), ZDim.end(),
              GravityTest.GetZCellSizes().begin()));

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

  //compare calculations for a single prims with known values
  BOOST_AUTO_TEST_CASE(box_gravity_calc_test)
    {
      //the gravity in the center of the cube should be 0
      //the density is 1 in all calculations here so do not have to multiply the term
      BOOST_CHECK(fabs(jiba::CalcGravBoxTerm(0.0, 0.0, 0.0, -10.0, -10.0, -10.0,
                  20.0, 20.0, 20.0)) < std::numeric_limits<double>::epsilon());
      //compare with the reported value of li and chouteau within a precision of 0.1%
      double topofcube = jiba::CalcGravBoxTerm(0.0, 0.0, 10.0, -10.0, -10.0,
          -10.0, 20.0, 20.0, 20.0);
      BOOST_CHECK_CLOSE(topofcube, -3.46426e-6, 0.1);
      double bottomofcube = jiba::CalcGravBoxTerm(0.0, 0.0, -10.0, -10.0, -10.0,
          -10.0, 20.0, 20.0, 20.0);
      //check symmetry of results
      BOOST_CHECK_CLOSE(topofcube, -bottomofcube,
          std::numeric_limits<float>::epsilon());
      double away1 = jiba::CalcGravBoxTerm(1e5, 10.0, 10.0, -10.0, -10.0, -10.0,
          20.0, 20.0, 20.0);
      BOOST_CHECK(away1 < std::numeric_limits<double>::epsilon());
      double away2 = jiba::CalcGravBoxTerm(100.0, 100.0, 10.0, -10.0, -10.0,
          -10.0, 20.0, 20.0, 20.0);
      BOOST_CHECK_CLOSE(away2, -1.873178e-9, 0.1);
      //check also for a very large cube
      BOOST_CHECK(fabs(jiba::CalcGravBoxTerm(0.0, 0.0, 0.0, -1e6, -1e6, -1e6,
                  2e6, 2e6, 2e6)) < std::numeric_limits<double>::epsilon());
    }

  //check that for a box the results are independent of the discretization
  BOOST_AUTO_TEST_CASE(model_gravity_boxcomp_test)
    {
      const double measx = 9.0;
      const double measy = 8.0;
      const double measz = -0.1;
      //again density is 1
      double boxtopofcube = jiba::CalcGravBoxTerm(measx, measy, measz, 0.0, 0.0,
          0.0, 20.0, 20.0, 20.0);
      jiba::GravimetryMatrix tensorbox = jiba::CalcTensorBoxTerm(measx, measy,
          measz, 0.0, 0.0, 0.0, 20.0, 20.0, 20.0);
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
      jiba::rvec DensityVector(ncells * ncells * ncells);
      for (size_t i = 0; i < ncells; ++i)
      for (size_t j = 0; j < ncells; ++j)
      for (size_t k = 0; k < ncells; ++k)
        {
          GravityTest.SetDensities()[i][j][k] = 1.0;
          DensityVector(i * (ncells * ncells) + j * ncells + k) = 1.0;
        }
      GravityTest.AddMeasurementPoint(measx, measy, measz);
      jiba::ThreeDGravityModel::tScalarMeasVec
      gravmeas(GravityTest.CalcGravity());
      jiba::ThreeDGravityModel::tTensorMeasVec tensmeas(
          GravityTest.CalcTensorGravity());
      double gridcube = gravmeas[0];

      GravityTest.WriteNetCDF("cube.nc");
      // Check that FTG calculation does not depend on discretization
      for (size_t i = 0; i < 2; ++i)
      for (size_t j = 0; j < 2; ++j)
        {
          BOOST_CHECK_CLOSE(tensmeas[0](i, j), tensorbox(i, j),
              std::numeric_limits<float>::epsilon());
        }

      //check some tensor properties
      BOOST_CHECK_CLOSE(tensmeas[0](0, 0) + tensmeas[0](1, 1),
          -tensmeas[0](2, 2), std::numeric_limits<float>::epsilon());
      BOOST_CHECK_CLOSE(tensmeas[0](1, 1) + tensmeas[0](2, 2),
          -tensmeas[0](0, 0), std::numeric_limits<float>::epsilon());
      BOOST_CHECK_CLOSE(tensmeas[0](0, 0) + tensmeas[0](2, 2),
          -tensmeas[0](1, 1), std::numeric_limits<float>::epsilon());
      BOOST_CHECK_CLOSE(boxtopofcube, gridcube,
          std::numeric_limits<float>::epsilon());
      //Check scalar sensitivity calculation
      jiba::rmat ScalarSensitivities(GravityTest.GetScalarSensitivities());
      jiba::rvec SensValues(boost::numeric::ublas::prec_prod(ScalarSensitivities,
              DensityVector));
      BOOST_CHECK_CLOSE(gridcube, SensValues(0),
          std::numeric_limits<float>::epsilon());
    }

  //check that the 1D background calculation gives correct results
  BOOST_AUTO_TEST_CASE(background_test)
    {
      jiba::ThreeDGravityModel GravityTest;
      double analytic = 2 * M_PI * jiba::Grav_const * (500.0 + 5.0 * 4500.0);
      std::vector<double> bg_dens, bg_thick;
      const size_t nmeas = 10;
      for (size_t i = 0; i < nmeas; ++i)
      GravityTest.AddMeasurementPoint(rand() % 50000 + 2e4, rand() % 50000
          + 2e4, 0.0);
      bg_dens += 1.0, 1.0, 5.0, 5.0;
      bg_thick += 200.0, 300.0, 3500.0, 1000.0;
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

  //compare the 3D solution with the analytic solution for a 1D layered structure
  BOOST_AUTO_TEST_CASE(model_gravity_1danaly_test)
    {
      jiba::ThreeDGravityModel GravityTest(true);
      //create a model with a random number of cells in each direction
      const size_t nhorcells = rand() % 50 + 10;
      const size_t nzcells = 10;
      GravityTest.SetXCellSizes().resize(boost::extents[nhorcells]);
      GravityTest.SetYCellSizes().resize(boost::extents[nhorcells]);
      GravityTest.SetZCellSizes().resize(boost::extents[nzcells]);
      jiba::rvec DensityVector(nhorcells * nhorcells * nzcells + 4); // 4 background layers
      for (size_t i = 0; i < nhorcells; ++i) // set the values of the inner cells

        {
          GravityTest.SetXCellSizes()[i] = rand() % 10000 + 1000;
          GravityTest.SetYCellSizes()[i] = rand() % 10000 + 1000;
        }
      for (size_t i = 0; i < nzcells; ++i)
        {
          GravityTest.SetZCellSizes()[i] = 500;
        }

      GravityTest.SetDensities().resize(
          boost::extents[nhorcells][nhorcells][nzcells]);
      for (size_t i = 0; i < nhorcells; ++i)
      for (size_t j = 0; j < nhorcells; ++j)
        {
          GravityTest.SetDensities()[i][j][0] = 1.0;
          DensityVector(i * nhorcells * nzcells + j * nzcells) = 1.0;
        }
      for (size_t i = 0; i < nhorcells; ++i)
      for (size_t j = 0; j < nhorcells; ++j)
      for (size_t k = 1; k < nzcells; ++k)
        {
          GravityTest.SetDensities()[i][j][k] = 5.0;
          DensityVector(i * nhorcells * nzcells + j * nzcells + k) = 5.0;
        }

      const size_t nmeas = 10;
      for (size_t i = 0; i < nmeas; ++i)
      GravityTest.AddMeasurementPoint(rand() % 50000 + 2e4, rand() % 50000
          + 2e4, 0.0);
      std::vector<double> bg_dens, bg_thick;
      bg_dens += 1.0, 1.0, 5.0, 5.0;
      bg_thick += 200.0, 300.0, 3500.0, 2000.0;

      DensityVector(nhorcells * nhorcells * nzcells) = 1.0;
      DensityVector(nhorcells * nhorcells * nzcells + 1) = 1.0;
      DensityVector(nhorcells * nhorcells * nzcells + 2) = 5.0;
      DensityVector(nhorcells * nhorcells * nzcells + 3) = 5.0;
      GravityTest.SetBackgroundDensities(bg_dens);
      GravityTest.SetBackgroundThicknesses(bg_thick);
      jiba::ThreeDGravityModel::tScalarMeasVec scalarmeas(
          GravityTest.CalcGravity());
      jiba::ThreeDGravityModel::tTensorMeasVec tensormeas(
          GravityTest.CalcTensorGravity());
      GravityTest.WriteNetCDF("layer.nc");
      double analytic = 2 * M_PI * jiba::Grav_const * (500.0 + 5.0 * 5500.0);
      jiba::rmat ScalarSensitivities(GravityTest.GetScalarSensitivities());
      jiba::rvec SensValues(prec_prod(ScalarSensitivities, DensityVector));
      //The second time the measurements will also be calculated from the sensitivities internally
      jiba::ThreeDGravityModel::tScalarMeasVec scalarmeas2(
          GravityTest.CalcGravity());
      for (size_t i = 0; i < nmeas; ++i)
        {
          BOOST_CHECK_CLOSE(analytic, scalarmeas[i], 0.01);
          BOOST_CHECK_CLOSE(scalarmeas[i], scalarmeas2[i], 0.001);
          BOOST_CHECK_CLOSE(analytic, SensValues[i], 0.01);
          BOOST_CHECK_CLOSE(scalarmeas[0], scalarmeas[i], 0.01);
          // The  tensorial elements should all be zero
          for (size_t i = 0; i < 2; ++i)
          for (size_t j = 0; j < 2; ++j)
            {
              BOOST_CHECK(fabs(tensormeas[i](i, j)) < std::numeric_limits<
                  double>::epsilon());
            }
        }
    }

  //check whether calculation by the sensitivity matrix yields the same results
  //for scalar calculations
  BOOST_AUTO_TEST_CASE(scalar_caching_test)
    {
      jiba::ThreeDGravityModel GravityTest(true, false);
      const size_t nmeas = 10;
      MakeRandomModel(GravityTest, nmeas);

      //Calculate twice, once with normal calculation, once cached
      //and record the times
      boost::posix_time::ptime startfirst =
      boost::posix_time::microsec_clock::local_time();
      jiba::ThreeDGravityModel::tScalarMeasVec scalarmeas1(
          GravityTest.CalcGravity());
      boost::posix_time::ptime endfirst =
      boost::posix_time::microsec_clock::local_time();
      jiba::ThreeDGravityModel::tScalarMeasVec scalarmeas2(
          GravityTest.CalcGravity());
      boost::posix_time::ptime endsecond =
      boost::posix_time::microsec_clock::local_time();
      //the second time it should be much faster
      BOOST_CHECK(endfirst - startfirst> (endsecond - endfirst));
      for (size_t i = 0; i < nmeas; ++i)
        {
          BOOST_CHECK_CLOSE(scalarmeas1[i], scalarmeas2[i], std::numeric_limits<
              float>::epsilon());
        }
      jiba::ScalarOMPGravityImp Implementation;
      jiba::rvec NewResult(jiba::MinMemGravityCalculator().Calculate(GravityTest,Implementation));
      for (size_t i = 0; i < nmeas; ++i)
        {
          BOOST_CHECK_CLOSE(scalarmeas1[i], NewResult[i], std::numeric_limits<
              float>::epsilon());
        }
    }

  //check some general properties of the FTG tensor that hold for any measurement above the surface
  BOOST_AUTO_TEST_CASE(random_tensor_test)
    {
      jiba::ThreeDGravityModel GravityTest;
      const size_t nmeas = 10;
      MakeRandomModel(GravityTest, nmeas);

      jiba::ThreeDGravityModel::tTensorMeasVec tensormeas(
          GravityTest.CalcTensorGravity());
      for (size_t i = 0; i < nmeas; ++i)
        {
          // check that tensor is traceless
          BOOST_CHECK_CLOSE(tensormeas[i](0, 0) + tensormeas[i](1, 1),
              -tensormeas[i](2, 2), std::numeric_limits<float>::epsilon());
          //check that tensor is symmetric
          BOOST_CHECK_CLOSE(tensormeas[i](0, 1), tensormeas[i](1, 0),
              std::numeric_limits<float>::epsilon());
          BOOST_CHECK_CLOSE(tensormeas[i](0, 2), tensormeas[i](2, 0),
              std::numeric_limits<float>::epsilon());
          BOOST_CHECK_CLOSE(tensormeas[i](1, 2), tensormeas[i](2, 1),
              std::numeric_limits<float>::epsilon());
        }
    }

  //compare with the ubc gravity forward code
  //BOOST_AUTO_TEST_CASE(ubc_test)
  //  {
  //    //create a random model
  //    jiba::ThreeDGravityModel GravityTest;
  //    const size_t nmeas = 10;
  //    MakeRandomModel(GravityTest, nmeas);
  //    //and compute the results for our code
  //    jiba::ThreeDGravityModel::tScalarMeasVec scalarmeas(
  //        GravityTest.CalcGravity());
  //    //set the mesh for the ubc code, they have the x-axis in east direction
  //    std::ofstream meshfile("testmesh");
  //    meshfile << GravityTest.GetYCellSizes().size() << " ";
  //    meshfile << GravityTest.GetXCellSizes().size() << " ";
  //    meshfile << GravityTest.GetZCellSizes().size() << "\n";
  //    meshfile << " 0 0 0 \n";
  //    //so we have to swap x and y in the grid coordinate system
  //    std::copy(GravityTest.GetYCellSizes().begin(),
  //        GravityTest.GetYCellSizes().end(), std::ostream_iterator<double>(
  //            meshfile, " "));
  //    meshfile << "\n";
  //    std::copy(GravityTest.GetXCellSizes().begin(),
  //        GravityTest.GetXCellSizes().end(), std::ostream_iterator<double>(
  //            meshfile, " "));
  //    meshfile << "\n";
  //    std::copy(GravityTest.GetZCellSizes().begin(),
  //        GravityTest.GetZCellSizes().end(), std::ostream_iterator<double>(
  //            meshfile, " "));
  //    meshfile << "\n" << std::flush;
  //    //write the observation points
  //    //again we have to switch x and y and z is positive up
  //    std::ofstream obsfile("testobs");
  //    obsfile << nmeas << "\n";
  //    for (size_t i = 0; i < nmeas; ++i)
  //      {
  //        obsfile << GravityTest.GetMeasPosY().at(i) << " "
  //            << GravityTest.GetMeasPosX().at(i) << " "
  //            << -GravityTest.GetMeasPosZ().at(i) << "\n";
  //      }
  //    obsfile << std::flush;
  //    //write out the densities, we can write them in the same order we store them
  //    std::ofstream densfile("testdens");
  //    std::copy(GravityTest.GetDensities().origin(),
  //        GravityTest.GetDensities().origin()
  //            + GravityTest.GetDensities().num_elements(),
  //        std::ostream_iterator<double>(densfile, "\n"));
  //    densfile << std::flush;
  //    //call the ubc code using wine
  //    system("wine gzfor3d.exe testmesh testobs testdens");
  //    //read in the results
  //    std::ifstream resultfile("gzfor3d.grv");
  //    //the first line contains some extra information we don't want
  //    char dummy[255];
  //    resultfile.getline(dummy, 255);
  //    //read in all numbers in the file into a vector
  //    //this includes coordinate information etc.
  //    std::vector<double> ubcresults;
  //    std::copy(std::istream_iterator<double>(resultfile),
  //        std::istream_iterator<double>(), std::back_inserter(ubcresults));
  //    for (size_t i = 0; i < nmeas; ++i)
  //      {
  //        //compare our results with the right component of the vector
  //        //the tolerance of 0.1% is OK considering truncation effect from writing to ascii
  //        BOOST_CHECK_CLOSE(scalarmeas[i] * 100000.0,
  //            ubcresults[(i + 1) * 4 - 1], 0.1);
  //      }
  //  }

  //compare the result from the tensor calculation with finite difference scalar calculations
  BOOST_AUTO_TEST_CASE(fd_tensor_test)
    {
      jiba::ThreeDGravityModel GravityTest;
      //we want to control the position of the measurements ourselves
      const size_t nmeas = 0;
      const double delta = 0.0001;
      //the value of delta in %
      const double precision = 0.05;
      MakeRandomModel(GravityTest, nmeas);
      //setup points for finite differencing in vertical and horizontal directions
      GravityTest.AddMeasurementPoint(50, 70, -5);
      GravityTest.AddMeasurementPoint(50, 70, -5 - delta);
      GravityTest.AddMeasurementPoint(50, 70, -5 + delta);
      GravityTest.AddMeasurementPoint(50, 70 - delta, -5);
      GravityTest.AddMeasurementPoint(50, 70 + delta, -5);
      GravityTest.AddMeasurementPoint(50 - delta, 70, -5);
      GravityTest.AddMeasurementPoint(50 + delta, 70, -5);
      //perform the calculations
      jiba::ThreeDGravityModel::tTensorMeasVec tensormeas(
          GravityTest.CalcTensorGravity());
      jiba::ThreeDGravityModel::tScalarMeasVec scalarmeas(
          GravityTest.CalcGravity());
      //extract the tensor components
      double uzz = tensormeas.at(1)(2, 2);
      double uzy = tensormeas.at(1)(2, 1);
      double uzx = tensormeas.at(1)(2, 0);
      //calculate the same components through finite differencing
      double fduzz = (scalarmeas.at(2) - scalarmeas.at(1)) / (2 * delta);
      double fduzy = (scalarmeas.at(4) - scalarmeas.at(3)) / (2 * delta);
      double fduzx = (scalarmeas.at(6) - scalarmeas.at(5)) / (2 * delta);
      //check that they match within the precision of delta
      BOOST_CHECK_CLOSE(uzz, fduzz, precision);
      BOOST_CHECK_CLOSE(uzy, fduzy, precision);
      BOOST_CHECK_CLOSE(uzx, fduzx, precision);
    }

  //check whether calculation by the sensitivity matrix yields the same results
  //for tensor calculations
  BOOST_AUTO_TEST_CASE(tensor_caching_test)
    {
      jiba::ThreeDGravityModel GravityTest(false, true);
      const size_t nmeas = 10;
      MakeRandomModel(GravityTest, nmeas);

      //Calculate twice, once with normal calculation, once cached
      boost::posix_time::ptime startfirst =
      boost::posix_time::microsec_clock::local_time();
      jiba::ThreeDGravityModel::tTensorMeasVec tensormeas1(
          GravityTest.CalcTensorGravity());
      boost::posix_time::ptime endfirst =
      boost::posix_time::microsec_clock::local_time();
      jiba::ThreeDGravityModel::tTensorMeasVec tensormeas2(
          GravityTest.CalcTensorGravity());
      boost::posix_time::ptime endsecond =
      boost::posix_time::microsec_clock::local_time();
      //the second time it should be much faster
      BOOST_CHECK(endfirst - startfirst> (endsecond - endfirst));
      //compare all tensor elements
      for (size_t i = 0; i < nmeas; ++i)
        {
          BOOST_CHECK_CLOSE(tensormeas1[i](0, 0), tensormeas2[i](0, 0),
              std::numeric_limits<float>::epsilon());
          BOOST_CHECK_CLOSE(tensormeas1[i](1, 0), tensormeas2[i](1, 0),
              std::numeric_limits<float>::epsilon());
          BOOST_CHECK_CLOSE(tensormeas1[i](2, 0), tensormeas2[i](2, 0),
              std::numeric_limits<float>::epsilon());
          BOOST_CHECK_CLOSE(tensormeas1[i](0, 1), tensormeas2[i](0, 1),
              std::numeric_limits<float>::epsilon());
          BOOST_CHECK_CLOSE(tensormeas1[i](1, 1), tensormeas2[i](1, 1),
              std::numeric_limits<float>::epsilon());
          BOOST_CHECK_CLOSE(tensormeas1[i](2, 1), tensormeas2[i](2, 1),
              std::numeric_limits<float>::epsilon());
          BOOST_CHECK_CLOSE(tensormeas1[i](0, 2), tensormeas2[i](0, 2),
              std::numeric_limits<float>::epsilon());
          BOOST_CHECK_CLOSE(tensormeas1[i](1, 2), tensormeas2[i](1, 2),
              std::numeric_limits<float>::epsilon());
          BOOST_CHECK_CLOSE(tensormeas1[i](2, 2), tensormeas2[i](2, 2),
              std::numeric_limits<float>::epsilon());
        }
    }

  //write a C++ vectorial quantity into a file so that R can understand the values
  template<typename VectorType>
  void WriteVectorToScript(std::ofstream &file, VectorType thevector,
      std::string Name)
    {
      file << Name << "<-c(";
      copy(thevector.begin(), thevector.end() - 1, std::ostream_iterator<double>(
              file, ","));
      file << thevector[thevector.size() - 1];
      file << ")\n";
    }

  //Make a random model without background layers and write the information
  //about the model into a file so that R can understand it
  void PrepareModelForR(jiba::ThreeDGravityModel &GravityTest,
      std::ofstream &scriptfile)
    {
      //first we setup our local model with some test measurements
      const size_t nmeas = 10;
      MakeRandomModel(GravityTest, nmeas);
      //we generate a script file for R that produces the same model
      //the build environment has copied these libraries in the right path
      scriptfile << "dyn.load(\"libmodelbase.so\") \n";
      scriptfile << "dyn.load(\"libgravity.so\") \n";
      scriptfile << "source(\"Gravity/GravForward.R\") \n";
      scriptfile << "alloc<-.C(\"AllocateModel\",as.integer(1),as.integer(1)) \n";
      //write the different model quantities into the script
      //the cell size coordinates
      WriteVectorToScript(scriptfile, GravityTest.GetXCellSizes(), "XSizes");
      WriteVectorToScript(scriptfile, GravityTest.GetYCellSizes(), "YSizes");
      WriteVectorToScript(scriptfile, GravityTest.GetZCellSizes(), "ZSizes");
      //and the measurement coordinates
      WriteVectorToScript(scriptfile, GravityTest.GetMeasPosX(), "XMeasPos");
      WriteVectorToScript(scriptfile, GravityTest.GetMeasPosY(), "YMeasPos");
      WriteVectorToScript(scriptfile, GravityTest.GetMeasPosZ(), "ZMeasPos");
      //copy the 3D density model into a vector
      jiba::rvec DensityVector(GravityTest.GetDensities().num_elements());
      copy(GravityTest.GetDensities().origin(),
          GravityTest.GetDensities().origin()
          + GravityTest.GetDensities().num_elements(), DensityVector.begin());
      //write the density vector into the script as well
      WriteVectorToScript(scriptfile, DensityVector, "Densities");
    }

  //test the scalar forward interface for R
  //this needs the current svn of boost test, as boost test in 1.35.0 has problems with the system call
  BOOST_AUTO_TEST_CASE(R_scalar_interface_test)
    {
      //create a 3D Gravity object
      jiba::ThreeDGravityModel GravityTest(false, false);
      //create a random model and write the information into the R-script file
      std::ofstream Rscript("scalar_test.R");
      PrepareModelForR(GravityTest, Rscript);
      //calculate our results
      jiba::ThreeDGravityModel::tScalarMeasVec scalarmeas(
          GravityTest.CalcGravity());
      //finish the R script
      //call the R interface function
      Rscript
      << " raw<-system.time(result<-gravforward(XSizes,YSizes,ZSizes,Densities,XMeasPos,YMeasPos,ZMeasPos))\n";
      Rscript
      << " cached<-system.time(result2<-gravforward(XSizes,YSizes,ZSizes,Densities,XMeasPos,YMeasPos,ZMeasPos))\n";
      Rscript << " sink(\"scalar_output\")\n";
      Rscript << " cat(result$GravAcceleration)\n";
      Rscript << " sink(\"r_timing\")\n";
      Rscript << " cat(raw)\n";
      Rscript << " cat(cached)\n";
      Rscript << " q()\n";
      Rscript << std::flush;
      //execute R with the script
      system("R --vanilla --slave -f scalar_test.R");
      //read in the output R has generated
      std::ifstream routput("scalar_output");
      std::vector<double> rvalues;
      //this is simply an ascii file with a bunch of numbers
      copy(std::istream_iterator<double>(routput), std::istream_iterator<double>(),
          back_inserter(rvalues));
      //first of all we should have values for each measurement
      BOOST_CHECK_EQUAL(rvalues.size(), scalarmeas.size());
      //and they should be equal, the tolerance is 0.01% as writing to the file truncates the numbers
      for (size_t i = 0; i < scalarmeas.size(); ++i)
        {
          BOOST_CHECK_CLOSE(scalarmeas.at(i), rvalues.at(i), 0.01);
        }

    }

  //test the tensor forward interface for R
  //this needs the current svn of boost test, as boost test in 1.35.0 has problems with the system call
  BOOST_AUTO_TEST_CASE(R_tensor_interface_test)
    {
      jiba::ThreeDGravityModel GravityTest(false, false);

      std::ofstream Rscript("tensor_test.R");
      PrepareModelForR(GravityTest, Rscript);
      jiba::ThreeDGravityModel::tTensorMeasVec tensormeas(
          GravityTest.CalcTensorGravity());
      Rscript
      << " result<-gravtensorforward(XSizes,YSizes,ZSizes,Densities,XMeasPos,YMeasPos,ZMeasPos)\n";
      Rscript << " sink(\"tensor_output\")\n";
      Rscript << " cat(result$GravAcceleration)\n";
      Rscript << " q()\n";
      Rscript << std::flush;
      //execute R with the script
      system("R --vanilla --slave -f tensor_test.R");
      //read in the output R has generated
      std::ifstream routput("tensor_output");
      std::vector<double> rvalues;
      //this is simply an ascii file with a bunch of numbers
      copy(std::istream_iterator<double>(routput), std::istream_iterator<double>(),
          back_inserter(rvalues));
      //first of all we should have values for each measurement
      //tensormeas is a vector of matrices while rvalues is just a vector where 9 consecutive elements
      //correspond to 1 FTG matrix
      BOOST_CHECK_EQUAL(rvalues.size(), tensormeas.size() * 9);
      //and they should be equal, the tolerance is 0.01% as writing to the file truncates the numbers
      for (size_t i = 0; i < tensormeas.size(); ++i)
        {
          BOOST_CHECK_CLOSE(tensormeas.at(i)(0, 0), rvalues.at(i * 9), 0.01);
          BOOST_CHECK_CLOSE(tensormeas.at(i)(0, 1), rvalues.at(i * 9 + 1), 0.01);
          BOOST_CHECK_CLOSE(tensormeas.at(i)(0, 2), rvalues.at(i * 9 + 2), 0.01);
          BOOST_CHECK_CLOSE(tensormeas.at(i)(1, 0), rvalues.at(i * 9 + 3), 0.01);
          BOOST_CHECK_CLOSE(tensormeas.at(i)(1, 1), rvalues.at(i * 9 + 4), 0.01);
          BOOST_CHECK_CLOSE(tensormeas.at(i)(1, 2), rvalues.at(i * 9 + 5), 0.01);
          BOOST_CHECK_CLOSE(tensormeas.at(i)(2, 0), rvalues.at(i * 9 + 6), 0.01);
          BOOST_CHECK_CLOSE(tensormeas.at(i)(2, 1), rvalues.at(i * 9 + 7), 0.01);
          BOOST_CHECK_CLOSE(tensormeas.at(i)(2, 2), rvalues.at(i * 9 + 8), 0.01);
        }

    }
  BOOST_AUTO_TEST_SUITE_END()
