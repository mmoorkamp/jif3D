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

  //Test the access function for the cell dimensions
  BOOST_AUTO_TEST_CASE(cell_dimension_access_function_test)
    {
      jiba::ThreeDGravityModel BaseTest;
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

  BOOST_AUTO_TEST_CASE(coordinate_calculation_test)
    {
      const size_t nelements = 5;
      jiba::ThreeDGravityModel BaseTest;
      jiba::ThreeDModelBase::t3DModelDim TestDim(boost::extents[nelements]);
      for (size_t i = 0; i < nelements; ++i)
        TestDim[i] = rand();
      BaseTest.SetXCellSizes().resize(boost::extents[nelements]);
      BaseTest.SetXCellSizes() = TestDim;
      BaseTest.SetYCellSizes().resize(boost::extents[nelements]);
      BaseTest.SetYCellSizes() = TestDim;
      BaseTest.SetZCellSizes().resize(boost::extents[nelements]);
      BaseTest.SetZCellSizes() = TestDim;
      jiba::ThreeDModelBase::t3DModelDim XCoordinates(BaseTest.GetXCoordinates());
      jiba::ThreeDModelBase::t3DModelDim YCoordinates(BaseTest.GetYCoordinates());
      jiba::ThreeDModelBase::t3DModelDim ZCoordinates(BaseTest.GetZCoordinates());
      BOOST_CHECK_EQUAL(XCoordinates[0],0.0);
      BOOST_CHECK_EQUAL(YCoordinates[0],0.0);
      BOOST_CHECK_EQUAL(ZCoordinates[0],0.0);
      jiba::ThreeDModelBase::t3DModelDim CompCoord(boost::extents[nelements]);
      std::partial_sum(TestDim.begin(), TestDim.end(), CompCoord.begin());
      BOOST_CHECK(std::equal(CompCoord.begin(), CompCoord.end()-1,
              BaseTest.GetXCoordinates().begin()+1));
      BOOST_CHECK(std::equal(CompCoord.begin(), CompCoord.end()-1,
              BaseTest.GetYCoordinates().begin()+1));
      BOOST_CHECK(std::equal(CompCoord.begin(), CompCoord.end()-1,
              BaseTest.GetZCoordinates().begin()+1));
    }

  BOOST_AUTO_TEST_CASE(find_indices_test)
    {
      const size_t nelements = 5;
      jiba::ThreeDGravityModel BaseTest;
      BaseTest.SetXCellSizes().resize(boost::extents[nelements]);
      BaseTest.SetYCellSizes().resize(boost::extents[nelements]);
      BaseTest.SetZCellSizes().resize(boost::extents[nelements]);
      for (size_t i = 0; i < nelements; ++i)
        {
          BaseTest.SetXCellSizes()[i] = 1.5;
          BaseTest.SetYCellSizes()[i] = 2.0;
          BaseTest.SetZCellSizes()[i] = 2.5;

        }
      boost::array<jiba::ThreeDModelBase::t3DModelData::index, 3> indices(
          BaseTest.FindAssociatedIndices(2.2, 2.2, 2.2));
      BOOST_CHECK_EQUAL(indices[0], 1);
      BOOST_CHECK_EQUAL(indices[1], 1);
      BOOST_CHECK_EQUAL(indices[2], 0);
    }


  //Check whether storing the model in a netcdf file and reading it in again creates the same model
  BOOST_AUTO_TEST_CASE(netcdf_read_write_test)
    {
      jiba::ThreeDGravityModel GravityTest;

      srand(time(NULL));

      jiba::ThreeDModelBase::t3DModelDim XDim = GenerateDimension(11);
      jiba::ThreeDModelBase::t3DModelDim YDim = GenerateDimension(23);
      jiba::ThreeDModelBase::t3DModelDim ZDim = GenerateDimension(14);
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
