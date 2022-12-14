//============================================================================
// Name        : test_ReadWriteX3D.cpp
// Author      : Feb 17, 2009
// Version     :
// Copyright   : 2009, mmoorkamp
//============================================================================

#define BOOST_TEST_MODULE ReadWriteX3D test
#define BOOST_TEST_MAIN ...
#include "../Global/Jif3DTesting.h"
#include <boost/test/included/unit_test.hpp>
#include <boost/test/tools/floating_point_comparison.hpp>
#include "ReadWriteX3D.h"
#include "../ModelBase/VTKTools.h"
#include "../ModelBase/NetCDFModelTools.h"
#include "../Global/Jif3DPlatformHelper.h"

BOOST_AUTO_TEST_SUITE( ReadWriteX3D_Suite )

    BOOST_AUTO_TEST_CASE (read_write_X3D_test)
      {
        //create a random number of cells and background layers
        const size_t xsize = rand() % 10 + 2;
        const size_t ysize = rand() % 10 + 2;
        const size_t zsize = rand() % 10 + 2;
        const size_t nbglayers = rand() % 10 + 2;
        std::vector<double> Depths(2);
        std::fill(Depths.begin(), Depths.end(), 12.3);
        jif3D::ThreeDModelBase::t3DModelDim XCellSizes(xsize), YCellSizes(ysize),
            ZCellSizes(zsize);
        jif3D::ThreeDModelBase::t3DModelData Data(boost::extents[xsize][ysize][zsize]);
        std::vector<double> bg_thicknesses(nbglayers), bg_conductivities(nbglayers);

        const double deltax = 11.2;
        const double deltay = 17.4;
        std::fill_n(XCellSizes.begin(), xsize, deltax);
        std::fill_n(YCellSizes.begin(), ysize, deltay);
        std::generate_n(ZCellSizes.begin(), zsize, jif3D::platform::drand48);
        std::generate_n(Data.origin(), xsize * ysize * zsize, jif3D::platform::drand48);
        std::generate_n(bg_thicknesses.begin(), nbglayers, jif3D::platform::drand48);
        std::generate_n(bg_conductivities.begin(), nbglayers, jif3D::platform::drand48);

        const std::string filename("x3d.in");
        jif3D::Write3DModelForX3D(filename, XCellSizes, YCellSizes, ZCellSizes, Depths,
            Data, bg_thicknesses, bg_conductivities);

        jif3D::ThreeDModelBase::t3DModelDim InXCellSizes, InYCellSizes, InZCellSizes;
        jif3D::ThreeDModelBase::t3DModelData InData;
        std::vector<double> Inbg_thicknesses, Inbg_conductivities;
        jif3D::Read3DModelFromX3D(filename, InXCellSizes, InYCellSizes, InZCellSizes,
            InData, Inbg_thicknesses, Inbg_conductivities);
        BOOST_CHECK(InXCellSizes.size() == xsize);
        BOOST_CHECK(InYCellSizes.size() == ysize);
        BOOST_CHECK(InZCellSizes.size() == zsize);
        BOOST_CHECK(Data.shape()[0] == xsize);
        BOOST_CHECK(Data.shape()[1] == ysize);
        BOOST_CHECK(Data.shape()[2] == zsize);
        BOOST_CHECK(Inbg_conductivities.size() == nbglayers);
        BOOST_CHECK(Inbg_thicknesses.size() == nbglayers);
        BOOST_CHECK(
            std::equal(XCellSizes.begin(), XCellSizes.end(), InXCellSizes.begin()));
        BOOST_CHECK(
            std::equal(YCellSizes.begin(), YCellSizes.end(), InYCellSizes.begin()));
        const double precision = 0.1;
        for (size_t i = 0; i < ZCellSizes.size(); ++i)
          {
            BOOST_CHECK_CLOSE(ZCellSizes[i], InZCellSizes[i], precision);
          }
        for (size_t i = 0; i < bg_thicknesses.size(); ++i)
          {
            BOOST_CHECK_CLOSE(bg_thicknesses[i], bg_thicknesses[i], precision);
          }
        for (size_t i = 0; i < bg_conductivities.size(); ++i)
          {
            BOOST_CHECK_CLOSE(bg_conductivities[i], Inbg_conductivities[i], precision);
          }
        for (size_t i = 0; i < Data.num_elements(); ++i)
          {
            BOOST_CHECK_CLOSE(*(Data.origin() + i), *(InData.origin() + i), precision);
          }

      }

    BOOST_AUTO_TEST_SUITE_END()
