//============================================================================
// Name        : test_readwritedata.cpp
// Author      : Jan 9, 2009
// Version     :
// Copyright   : 2009, mmoorkamp
//============================================================================

#define BOOST_TEST_MODULE ReadWriteGravityData test
#define BOOST_TEST_MAIN ...
#include "ReadWriteGravityData.h"
#include <boost/test/included/unit_test.hpp>
#include <boost/test/included/unit_test_framework.hpp>

#include <boost/test/floating_point_comparison.hpp>
#include <stdlib.h>

BOOST_AUTO_TEST_SUITE( GravityReadWrite_Test_Suite )

BOOST_AUTO_TEST_CASE(readwrite_scalar_test)
      {
        const size_t ndata = 20;
        jiba::rvec Data(ndata);
        jiba::ThreeDGravityModel::tMeasPosVec PosX(ndata, 0.0),
            PosY(ndata, 0.0), PosZ(ndata, 0.0);
        std::generate_n(Data.begin(), ndata, rand);
        std::generate_n(PosX.begin(), ndata, rand);
        std::generate_n(PosY.begin(), ndata, rand);
        std::generate_n(PosZ.begin(), ndata, rand);
        std::string filename("scalar.nc");
        jiba::SaveScalarGravityMeasurements(filename, Data, PosX, PosY, PosZ);
        jiba::rvec NewData(ndata);
        jiba::ThreeDGravityModel::tMeasPosVec NewPosX(ndata, 0.0), NewPosY(
            ndata, 0.0), NewPosZ(ndata, 0.0);
        jiba::ReadScalarGravityMeasurements(filename, NewData, NewPosX,
            NewPosY, NewPosZ);
        BOOST_CHECK(std::equal(Data.begin(), Data.end(),NewData.begin()));
        BOOST_CHECK(std::equal(PosX.begin(), PosX.end(),NewPosX.begin()));
        BOOST_CHECK(std::equal(PosY.begin(), PosY.end(),NewPosY.begin()));
        BOOST_CHECK(std::equal(PosZ.begin(), PosZ.end(),NewPosZ.begin()));
        BOOST_CHECK(jiba::IdentifyGravityDatafileType(filename) == jiba::scalar);
      }

BOOST_AUTO_TEST_CASE(readwrite_tensor_test)
      {
        const size_t nmeas  =20;
        const size_t ndata = nmeas*9;
        jiba::rvec Data(ndata);
        jiba::ThreeDGravityModel::tMeasPosVec PosX(nmeas, 0.0),
            PosY(nmeas, 0.0), PosZ(nmeas, 0.0);
        std::generate_n(Data.begin(), ndata, rand);
        std::generate_n(PosX.begin(), nmeas, rand);
        std::generate_n(PosY.begin(), nmeas, rand);
        std::generate_n(PosZ.begin(), nmeas, rand);
        std::string filename("tensor.nc");
        jiba::SaveTensorGravityMeasurements(filename, Data, PosX, PosY, PosZ);
        jiba::rvec NewData(ndata);
        jiba::ThreeDGravityModel::tMeasPosVec NewPosX(nmeas, 0.0), NewPosY(
            nmeas, 0.0), NewPosZ(nmeas, 0.0);
        jiba::ReadTensorGravityMeasurements(filename, NewData, NewPosX,
            NewPosY, NewPosZ);
        BOOST_CHECK(std::equal(Data.begin(), Data.end(),NewData.begin()));
        BOOST_CHECK(std::equal(PosX.begin(), PosX.end(),NewPosX.begin()));
        BOOST_CHECK(std::equal(PosY.begin(), PosY.end(),NewPosY.begin()));
        BOOST_CHECK(std::equal(PosZ.begin(), PosZ.end(),NewPosZ.begin()));
        BOOST_CHECK(jiba::IdentifyGravityDatafileType(filename) == jiba::ftg);
      }
    BOOST_AUTO_TEST_SUITE_END()
