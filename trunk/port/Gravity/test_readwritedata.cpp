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


//check whether writing data to a file and reading it in again
//gives the same data back for the scalar case
BOOST_AUTO_TEST_CASE(readwrite_scalar_test)
      {
        const size_t ndata = 20;
        jif3D::rvec Data(ndata), Error(ndata);
        jif3D::ThreeDGravityModel::tMeasPosVec PosX(ndata, 0.0),
            PosY(ndata, 0.0), PosZ(ndata, 0.0);
        std::generate_n(Data.begin(), ndata, rand);
        std::generate_n(Error.begin(), ndata, rand);
        std::generate_n(PosX.begin(), ndata, rand);
        std::generate_n(PosY.begin(), ndata, rand);
        std::generate_n(PosZ.begin(), ndata, rand);
        std::string filename("scalar.nc");
        jif3D::SaveScalarGravityMeasurements(filename, Data, PosX, PosY, PosZ, Error);
        jif3D::rvec NewData(ndata), NewError(ndata);
        jif3D::ThreeDGravityModel::tMeasPosVec NewPosX(ndata, 0.0), NewPosY(
            ndata, 0.0), NewPosZ(ndata, 0.0);
        jif3D::ReadScalarGravityMeasurements(filename, NewData, NewPosX,
            NewPosY, NewPosZ, NewError);
        BOOST_CHECK(std::equal(Data.begin(), Data.end(),NewData.begin()));
        BOOST_CHECK(std::equal(Error.begin(), Error.end(),Error.begin()));
        BOOST_CHECK(std::equal(PosX.begin(), PosX.end(),NewPosX.begin()));
        BOOST_CHECK(std::equal(PosY.begin(), PosY.end(),NewPosY.begin()));
        BOOST_CHECK(std::equal(PosZ.begin(), PosZ.end(),NewPosZ.begin()));
        BOOST_CHECK(jif3D::IdentifyGravityDatafileType(filename) == jif3D::scalar);
      }

//check whether writing data to a file and reading it in again
//gives the same data back for the tensor case
BOOST_AUTO_TEST_CASE(readwrite_tensor_test)
      {
        const size_t nmeas  =20;
        const size_t ndata = nmeas*9;
        jif3D::rvec Data(ndata), Error(ndata);
        jif3D::ThreeDGravityModel::tMeasPosVec PosX(nmeas, 0.0),
            PosY(nmeas, 0.0), PosZ(nmeas, 0.0);
        std::generate_n(Data.begin(), ndata, rand);
        std::generate_n(Error.begin(), ndata, rand);
        std::generate_n(PosX.begin(), nmeas, rand);
        std::generate_n(PosY.begin(), nmeas, rand);
        std::generate_n(PosZ.begin(), nmeas, rand);
        std::string filename("tensor.nc");
        jif3D::SaveTensorGravityMeasurements(filename, Data, PosX, PosY, PosZ, Error);
        jif3D::rvec NewData(ndata), NewError(ndata);
        jif3D::ThreeDGravityModel::tMeasPosVec NewPosX(nmeas, 0.0), NewPosY(
            nmeas, 0.0), NewPosZ(nmeas, 0.0);
        jif3D::ReadTensorGravityMeasurements(filename, NewData, NewPosX,
            NewPosY, NewPosZ, NewError);
        BOOST_CHECK(std::equal(Data.begin(), Data.end(),NewData.begin()));
        BOOST_CHECK(std::equal(Error.begin(), Error.end(),NewError.begin()));
        BOOST_CHECK(std::equal(PosX.begin(), PosX.end(),NewPosX.begin()));
        BOOST_CHECK(std::equal(PosY.begin(), PosY.end(),NewPosY.begin()));
        BOOST_CHECK(std::equal(PosZ.begin(), PosZ.end(),NewPosZ.begin()));
        BOOST_CHECK(jif3D::IdentifyGravityDatafileType(filename) == jif3D::ftg);
      }
    BOOST_AUTO_TEST_SUITE_END()
