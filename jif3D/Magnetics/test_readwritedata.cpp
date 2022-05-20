//============================================================================
// Name        : test_readwritedata.cpp
// Author      : Jan 9, 2009
// Version     :
// Copyright   : 2009, mmoorkamp
//============================================================================

#define BOOST_TEST_MODULE ReadWriteGravityData test
#define BOOST_TEST_MAIN ...
#include "../Global/Jif3DTesting.h"
#include "TotalFieldMagneticData.h"
#include "ThreeComponentMagneticData.h"
#include <boost/test/included/unit_test.hpp>
#include <boost/test/tools/floating_point_comparison.hpp>

BOOST_AUTO_TEST_SUITE (MagneticsReadWrite_Test_Suite)

//check whether writing data to a file and reading it in again
//gives the same data back for the scalar case
BOOST_AUTO_TEST_CASE(readwrite_total_test)
  {
    const size_t ndata = rand() % 200;
    std::vector<double> Data(ndata), Error(ndata);
    std::vector<double> PosX(ndata, 0.0), PosY(ndata, 0.0), PosZ(ndata, 0.0);
    std::generate_n(Data.begin(), ndata, rand);
    std::generate_n(Error.begin(), ndata, rand);
    std::generate_n(PosX.begin(), ndata, rand);
    std::generate_n(PosY.begin(), ndata, rand);
    std::generate_n(PosZ.begin(), ndata, rand);
    jif3D::TotalFieldMagneticData TotalData;
    TotalData.SetMeasurementPoints(PosX,PosY,PosZ);
    TotalData.SetDataAndErrors(Data,Error);
    std::string filename("total.nc");
    TotalData.WriteNetCDF(filename);
    jif3D::TotalFieldMagneticData CompData;
    CompData.ReadNetCDF("total.nc");

    BOOST_TEST(TotalData.GetData() == CompData.GetData(),
        boost::test_tools::per_element());
    BOOST_TEST(TotalData.GetErrors() == CompData.GetErrors(),
        boost::test_tools::per_element());
    BOOST_TEST(TotalData.GetMeasPosX() == CompData.GetMeasPosX(),
        boost::test_tools::per_element());
    BOOST_TEST(TotalData.GetMeasPosY() == CompData.GetMeasPosY(),
        boost::test_tools::per_element());
    BOOST_TEST(TotalData.GetMeasPosZ() == CompData.GetMeasPosZ(),
        boost::test_tools::per_element());
  }

//check whether writing data to a file and reading it in again
//gives the same data back for the tensor case
BOOST_AUTO_TEST_CASE(readwrite_threecomp_test)
  {
    const size_t ndata = rand() % 200;
       std::vector<double> Data(3 * ndata), Error(3 * ndata);
       std::vector<double> PosX(ndata, 0.0), PosY(ndata, 0.0), PosZ(ndata, 0.0);
       std::generate_n(Data.begin(), ndata, rand);
       std::generate_n(Error.begin(), ndata, rand);
       std::generate_n(PosX.begin(), ndata, rand);
       std::generate_n(PosY.begin(), ndata, rand);
       std::generate_n(PosZ.begin(), ndata, rand);
       jif3D::ThreeComponentMagneticData ThreecData;
       ThreecData.SetMeasurementPoints(PosX,PosY,PosZ);
       ThreecData.SetDataAndErrors(Data,Error);
       std::string filename("threec.nc");
       ThreecData.WriteNetCDF(filename);
       jif3D::ThreeComponentMagneticData CompData;
       CompData.ReadNetCDF("threec.nc");

       BOOST_TEST(ThreecData.GetData() == CompData.GetData(),
           boost::test_tools::per_element());
       BOOST_TEST(ThreecData.GetErrors() == CompData.GetErrors(),
           boost::test_tools::per_element());
       BOOST_TEST(ThreecData.GetMeasPosX() == CompData.GetMeasPosX(),
           boost::test_tools::per_element());
       BOOST_TEST(ThreecData.GetMeasPosY() == CompData.GetMeasPosY(),
           boost::test_tools::per_element());
       BOOST_TEST(ThreecData.GetMeasPosZ() == CompData.GetMeasPosZ(),
           boost::test_tools::per_element());
  }
BOOST_AUTO_TEST_SUITE_END()
