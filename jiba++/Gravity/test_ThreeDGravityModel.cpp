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
#include <boost/test/floating_point_comparison.hpp>
#include <boost/assign/std/vector.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/numeric/conversion/cast.hpp>

#include <stdlib.h>
#include <time.h>
#include <fstream>
#include "test_common.h"
#include "ThreeDGravityModel.h"

BOOST_AUTO_TEST_SUITE( GravityModel_Test_Suite )

//Test the default state of the object
BOOST_AUTO_TEST_CASE  (constructors_test)
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

      const size_t nbglayers = rand() % 50;
      std::vector<double> bg_thick, bg_dens;
      std::generate_n(back_inserter(bg_thick),nbglayers,rand);
      std::generate_n(back_inserter(bg_dens),nbglayers,rand);
      GravityTest.SetBackgroundDensities(bg_dens);
      GravityTest.SetBackgroundThicknesses(bg_thick);
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
      BOOST_CHECK(std::equal(bg_thick.begin(),bg_thick.end(),NetCDFReadTest.GetBackgroundThicknesses().begin()));
      BOOST_CHECK(std::equal(bg_dens.begin(),bg_dens.end(),NetCDFReadTest.GetBackgroundDensities().begin()));
      TestData[0][0][0] += 1.0;
      BOOST_CHECK(!std::equal(TestData.begin(), TestData.end(),
              NetCDFReadTest.GetDensities().begin()));
    }

  BOOST_AUTO_TEST_CASE(ascii_measpos_test)
    {
      std::ofstream outfile;
      const size_t nmeas = 20;
      std::string filename = "meas.test";
      outfile.open(filename.c_str());
      std::vector<double> posx,posy,posz;
      for (size_t i = 0; i < nmeas; ++i)
        {
          posx.push_back(rand()%1000);
          posy.push_back(rand()%1000);
          posz.push_back(rand()%1000);
          outfile << posx.at(i) << " " << posy.at(i) << " " << posz.at(i) << "\n";
        }
      outfile.close();
      jiba::ThreeDGravityModel TestModel;
      TestModel.ReadMeasPosAscii(filename);
      for (size_t i = 0; i < nmeas; ++i)
        {
          BOOST_CHECK_EQUAL(posx.at(i),TestModel.GetMeasPosX()[i]);
          BOOST_CHECK_EQUAL(posy.at(i),TestModel.GetMeasPosY()[i]);
          BOOST_CHECK_EQUAL(posz.at(i),TestModel.GetMeasPosZ()[i]);
        }
    }

  BOOST_AUTO_TEST_CASE(netcdf_measpos_test)
    {
      jiba::ThreeDGravityModel GravityTest;
      const size_t nmeas = 10;
      MakeRandomModel(GravityTest, nmeas);
    }

  BOOST_AUTO_TEST_SUITE_END()
