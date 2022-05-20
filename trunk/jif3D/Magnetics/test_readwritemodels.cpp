//============================================================================
// Name        : test_readwritedata.cpp
// Author      : Jan 9, 2009
// Version     :
// Copyright   : 2009, mmoorkamp
//============================================================================

#define BOOST_TEST_MODULE ReadWriteGravityData test
#define BOOST_TEST_MAIN ...
#include "../Global/Jif3DTesting.h"
#include "ThreeDMagnetizationModel.h"
#include "ThreeDSusceptibilityModel.h"

#include <boost/test/included/unit_test.hpp>
#include <boost/test/tools/floating_point_comparison.hpp>

BOOST_AUTO_TEST_SUITE (MagneticModelReadWrite_Test_Suite)

//check whether writing data to a file and reading it in again
//gives the same data back for the scalar case
    BOOST_AUTO_TEST_CASE(readwrite_susceptibility_test)
      {
        const size_t xsize = rand() % 10 + 2;
        const size_t ysize = rand() % 10 + 2;
        const size_t zsize = rand() % 10 + 2;
        const size_t nbglayers = rand() % 10 + 2;
        const size_t nmodel = xsize * ysize * zsize;
        jif3D::ThreeDModelBase::t3DModelDim XCellSizes(xsize), YCellSizes(ysize),
            ZCellSizes(zsize);
        jif3D::ThreeDModelBase::t3DModelData Data(boost::extents[xsize][ysize][zsize]);
        std::vector<double> bg_thicknesses(nbglayers), bg_susceptibilities(nbglayers);
        std::generate(XCellSizes.begin(), XCellSizes.end(), drand48);
        std::generate(YCellSizes.begin(), YCellSizes.end(), drand48);
        std::generate(ZCellSizes.begin(), ZCellSizes.end(), drand48);
        std::generate(Data.origin(), Data.origin() + nmodel, drand48);
        std::generate(bg_thicknesses.begin(), bg_thicknesses.end(), drand48);
        std::generate(bg_susceptibilities.begin(), bg_susceptibilities.end(), drand48);
        jif3D::ThreeDSusceptibilityModel Model;
        Model.SetMeshSize(xsize, ysize, zsize);
        Model.SetBackgroundSusceptibilities(bg_susceptibilities);
        Model.SetBackgroundThicknesses(bg_thicknesses);
        Model.SetData() = Data;
        Model.SetXCellSizes(XCellSizes);
        Model.SetYCellSizes(YCellSizes);
        Model.SetZCellSizes(ZCellSizes);

        Model.WriteNetCDF("sus.nc");
        jif3D::ThreeDSusceptibilityModel CompModel;
        CompModel.ReadNetCDF("sus.nc");
        BOOST_TEST(Model.GetData() == CompModel.GetData(),
            boost::test_tools::per_element());
        BOOST_TEST(Model.GetXCellSizes() == CompModel.GetXCellSizes(),
            boost::test_tools::per_element());
        BOOST_TEST(Model.GetYCellSizes() == CompModel.GetYCellSizes(),
            boost::test_tools::per_element());
        BOOST_TEST(Model.GetZCellSizes() == CompModel.GetZCellSizes(),
            boost::test_tools::per_element());
        BOOST_TEST(
            Model.GetBackgroundThicknesses() == CompModel.GetBackgroundThicknesses(),
            boost::test_tools::per_element());
        BOOST_TEST(
            Model.GetBackgroundSusceptibilities() == CompModel.GetBackgroundSusceptibilities(),
            boost::test_tools::per_element());
      }

//check whether writing data to a file and reading it in again
//gives the same data back for the tensor case
    BOOST_AUTO_TEST_CASE(readwrite_magnetization_test)
      {
        const size_t xsize = rand() % 10 + 2;
        const size_t ysize = rand() % 10 + 2;
        const size_t zsize = rand() % 10 + 2;
        const size_t nbglayers = rand() % 10 + 2;
        const size_t nmodel = xsize * ysize * zsize;
        jif3D::ThreeDModelBase::t3DModelDim XCellSizes(xsize), YCellSizes(ysize),
            ZCellSizes(zsize);
        jif3D::ThreeDModelBase::t3DModelData MagX(boost::extents[xsize][ysize][zsize]);
        jif3D::ThreeDModelBase::t3DModelData MagY(boost::extents[xsize][ysize][zsize]);
        jif3D::ThreeDModelBase::t3DModelData MagZ(boost::extents[xsize][ysize][zsize]);
        std::vector<double> bg_thicknesses(nbglayers), bg_susceptibilities(nbglayers);
        std::generate(XCellSizes.begin(), XCellSizes.end(), drand48);
        std::generate(YCellSizes.begin(), YCellSizes.end(), drand48);
        std::generate(ZCellSizes.begin(), ZCellSizes.end(), drand48);
        std::generate(MagX.origin(), MagX.origin() + nmodel, drand48);
        std::generate(MagY.origin(), MagY.origin() + nmodel, drand48);
        std::generate(MagZ.origin(), MagZ.origin() + nmodel, drand48);

        std::generate(bg_thicknesses.begin(), bg_thicknesses.end(), drand48);
        std::generate(bg_susceptibilities.begin(), bg_susceptibilities.end(), drand48);
        jif3D::ThreeDMagnetizationModel Model;
        Model.SetMeshSize(xsize, ysize, zsize);
        Model.SetMagnetization_Y().resize(boost::extents[xsize][ysize][zsize]);
        Model.SetMagnetization_Z().resize(boost::extents[xsize][ysize][zsize]);
        Model.SetMagnetization_X() = MagX;
        Model.SetMagnetization_Y() = MagY;
        Model.SetMagnetization_Z() = MagZ;

        Model.SetXCellSizes(XCellSizes);
        Model.SetYCellSizes(YCellSizes);
        Model.SetZCellSizes(ZCellSizes);

        Model.WriteNetCDF("mag.nc");
        jif3D::ThreeDMagnetizationModel CompModel;
        CompModel.ReadNetCDF("mag.nc");
        BOOST_TEST(Model.GetMagnetization_X() == CompModel.GetMagnetization_X(),
            boost::test_tools::per_element());
        BOOST_TEST(Model.GetMagnetization_Y() == CompModel.GetMagnetization_Y(),
            boost::test_tools::per_element());
        BOOST_TEST(Model.GetMagnetization_Z() == CompModel.GetMagnetization_Z(),
            boost::test_tools::per_element());
        BOOST_TEST(Model.GetXCellSizes() == CompModel.GetXCellSizes(),
            boost::test_tools::per_element());
        BOOST_TEST(Model.GetYCellSizes() == CompModel.GetYCellSizes(),
            boost::test_tools::per_element());
        BOOST_TEST(Model.GetZCellSizes() == CompModel.GetZCellSizes(),
            boost::test_tools::per_element());

      }
    BOOST_AUTO_TEST_SUITE_END()
