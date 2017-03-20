//============================================================================
// Name        : test_ReadWriteTitanData.cpp
// Author      : Feb 17, 2009
// Version     :
// Copyright   : 2009, mmoorkamp
//============================================================================

#define BOOST_TEST_MODULE ReadWriteTitanData test
#define BOOST_TEST_MAIN ...
#include "../Global/Jif3DTesting.h"
#include <boost/test/included/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

#include "../Global/Jif3DPlatformHelper.h"
#include "ReadWriteTitanData.h"

BOOST_AUTO_TEST_SUITE( ReadWriteTitanData_Suite )

    BOOST_AUTO_TEST_CASE (read_write_netcdf_test)
      {
        const size_t nfreq = 5;
        const size_t nmeas = 10;
        const size_t nstat = 7;
        const std::string filename("titan.nc");
        std::vector<double> Frequencies(nfreq);
        std::vector<double> XCoord(nmeas), YCoord(nmeas), ZCoord(nmeas);
        std::vector<int> ExIndices(nfreq * nstat), EyIndices(nfreq * nstat), HIndices(
            nfreq * nstat);
        std::vector<double> C(nstat * 4);
        const size_t ndata = nfreq * nstat * 8;
        jif3D::rvec Impedances(ndata), Error(ndata);

        std::generate_n(Frequencies.begin(), nfreq, jif3D::platform::drand48);
        std::sort(Frequencies.begin(), Frequencies.end());
        std::generate_n(XCoord.begin(), nmeas, jif3D::platform::drand48);
        std::generate_n(YCoord.begin(), nmeas, jif3D::platform::drand48);
        std::generate_n(ZCoord.begin(), nmeas, jif3D::platform::drand48);

        std::generate_n(ExIndices.begin(), nfreq * nstat, std::rand);
        std::generate_n(EyIndices.begin(), nfreq * nstat, std::rand);
        std::generate_n(HIndices.begin(), nfreq * nstat, std::rand);

        std::generate_n(Impedances.begin(), ndata, jif3D::platform::drand48);
        std::generate_n(Error.begin(), ndata, jif3D::platform::drand48);
        std::generate_n(C.begin(), nstat * 4, jif3D::platform::drand48);
        for (size_t i = 1; i < ndata; ++i)
          {
            Error(i) = Error(i - 1);
          }
        jif3D::WriteTitanDataToNetCDF(filename, Frequencies, XCoord, YCoord, ZCoord,
            ExIndices, EyIndices, HIndices, Impedances, Error, C);

        std::vector<double> ReadFrequencies;
        std::vector<double> ReadXCoord, ReadYCoord, ReadZCoord, ReadC;
        std::vector<int> ReadExIndices, ReadEyIndices, ReadHIndices;
        jif3D::rvec ReadImpedances, ReadError;
        jif3D::ReadTitanDataFromNetCDF(filename, ReadFrequencies, ReadXCoord, ReadYCoord,
            ReadZCoord, ReadExIndices, ReadEyIndices, ReadHIndices, ReadImpedances,
            ReadError, ReadC);
        for (size_t i = 0; i < nfreq; ++i)
          {
            BOOST_CHECK_CLOSE(Frequencies[i], ReadFrequencies[i], 0.001);
          }
        for (size_t i = 0; i < nstat; ++i)
          {
            BOOST_CHECK_CLOSE(XCoord[i], ReadXCoord[i], 0.001);
            BOOST_CHECK_CLOSE(YCoord[i], ReadYCoord[i], 0.001);
            BOOST_CHECK_CLOSE(ZCoord[i], ReadZCoord[i], 0.001);
          }

        for (size_t i = 0; i < nfreq * nstat; ++i)
          {
            BOOST_CHECK_EQUAL(ExIndices[i], ReadExIndices[i]);
            BOOST_CHECK_EQUAL(EyIndices[i], ReadEyIndices[i]);
            BOOST_CHECK_EQUAL(HIndices[i], ReadHIndices[i]);
          }

        for (size_t i = 0; i < ndata; ++i)
          {
            BOOST_CHECK_CLOSE(Impedances(i), ReadImpedances(i), 0.001);
            BOOST_CHECK_CLOSE(Error(i), ReadError(i), 0.001);
          }

        for (size_t i = 0; i < nstat * 4; ++i)
          {
            BOOST_CHECK_CLOSE(C[i], ReadC[i], 0.001);
          }
      }

    BOOST_AUTO_TEST_SUITE_END()
