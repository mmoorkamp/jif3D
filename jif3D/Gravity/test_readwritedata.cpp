//============================================================================
// Name        : test_readwritedata.cpp
// Author      : Jan 9, 2009
// Version     :
// Copyright   : 2009, mmoorkamp
//============================================================================

#define BOOST_TEST_MODULE ReadWriteGravityData test
#define BOOST_TEST_MAIN ...
#include "../Global/Jif3DTesting.h"
#include "ReadWriteGravityData.h"
#include "ScalarGravityData.h"
#include "TensorGravityData.h"
#include <boost/test/included/unit_test.hpp>
#include <boost/test/tools/floating_point_comparison.hpp>

BOOST_AUTO_TEST_SUITE (GravityReadWrite_Test_Suite)

//check whether writing data to a file and reading it in again
//gives the same data back for the scalar case
    BOOST_AUTO_TEST_CASE(readwrite_scalar_test)
      {
        const size_t ndata = 20;
        std::vector<double> Data(ndata), Error(ndata);
        std::vector<double> PosX(ndata, 0.0), PosY(ndata, 0.0), PosZ(ndata, 0.0);
        std::generate_n(Data.begin(), ndata, rand);
        std::generate_n(Error.begin(), ndata, rand);
        std::generate_n(PosX.begin(), ndata, rand);
        std::generate_n(PosY.begin(), ndata, rand);
        std::generate_n(PosZ.begin(), ndata, rand);
        jif3D::ScalarGravityData ScalData;
        ScalData.SetMeasurementPoints(PosX, PosY, PosZ);
        ScalData.SetDataAndErrors(Data, Error);
        std::string filename("scalar.nc");
        ScalData.WriteNetCDF(filename);
        jif3D::ScalarGravityData CompData;
        CompData.ReadNetCDF("scalar.nc");

        BOOST_TEST(ScalData.GetData() == CompData.GetData(),
            boost::test_tools::per_element());
        BOOST_TEST(ScalData.GetErrors() == CompData.GetErrors(),
            boost::test_tools::per_element());
        BOOST_TEST(ScalData.GetMeasPosX() == CompData.GetMeasPosX(),
            boost::test_tools::per_element());
        BOOST_TEST(ScalData.GetMeasPosY() == CompData.GetMeasPosY(),
            boost::test_tools::per_element());
        BOOST_TEST(ScalData.GetMeasPosZ() == CompData.GetMeasPosZ(),
            boost::test_tools::per_element());

        BOOST_CHECK(jif3D::IdentifyGravityDatafileType(filename) == jif3D::scalar);
      }

//check whether writing data to a file and reading it in again
//gives the same data back for the tensor case
    BOOST_AUTO_TEST_CASE(readwrite_tensor_test)
      {
        const size_t nmeas = 20;
        const size_t ndata = nmeas * 9;
        std::vector<double> Data(ndata), Error(ndata);
        std::vector<double> PosX(nmeas, 0.0), PosY(nmeas, 0.0), PosZ(nmeas, 0.0);
        std::generate_n(Data.begin(), ndata, rand);
        std::generate_n(Error.begin(), ndata, rand);
        std::generate_n(PosX.begin(), nmeas, rand);
        std::generate_n(PosY.begin(), nmeas, rand);
        std::generate_n(PosZ.begin(), nmeas, rand);
        jif3D::TensorGravityData TensData;
        TensData.SetMeasurementPoints(PosX, PosY, PosZ);
        TensData.SetDataAndErrors(Data, Error);
        std::string filename("tensor.nc");
        TensData.WriteNetCDF(filename);
        jif3D::TensorGravityData CompData;
        CompData.ReadNetCDF(filename);
        BOOST_TEST(TensData.GetData() == CompData.GetData(),
            boost::test_tools::per_element());
        BOOST_TEST(TensData.GetErrors() == CompData.GetErrors(),
            boost::test_tools::per_element());
        BOOST_TEST(TensData.GetMeasPosX() == CompData.GetMeasPosX(),
            boost::test_tools::per_element());
        BOOST_TEST(TensData.GetMeasPosY() == CompData.GetMeasPosY(),
            boost::test_tools::per_element());
        BOOST_TEST(TensData.GetMeasPosZ() == CompData.GetMeasPosZ(),
            boost::test_tools::per_element());
        BOOST_CHECK(jif3D::IdentifyGravityDatafileType(filename) == jif3D::ftg);
      }
    BOOST_AUTO_TEST_SUITE_END()
