//============================================================================
// Name        : test_ReadWriteX3D.cpp
// Author      : Feb 17, 2009
// Version     :
// Copyright   : 2009, mmoorkamp
//============================================================================

#define BOOST_TEST_MODULE ReadWriteImpedances test
#define BOOST_TEST_MAIN ...
#include "../Global/Jif3DTesting.h"
#include <boost/test/included/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

#include "../Global/Jif3DPlatformHelper.h"
#include "ReadWriteImpedances.h"
#include "X3DModel.h"
BOOST_AUTO_TEST_SUITE( ReadWriteImpedances_Suite )

    void GenerateData(std::vector<double> &Frequencies, std::vector<double> &XCoord,
        std::vector<double> &YCoord, std::vector<double> &ZCoord, jif3D::rvec &Impedances,
        jif3D::rvec &Error, std::vector<double> &C, size_t nstat = 7)
      {
        const size_t nfreq = 5;
        Frequencies.resize(nfreq);
        XCoord.resize(nstat);
        YCoord.resize(nstat);
        ZCoord.resize(nstat);
        const size_t ndata = nfreq * nstat * 8;
        Impedances.resize(ndata), Error.resize(ndata);
        C.resize(nstat * 4);
        std::generate_n(Frequencies.begin(), nfreq, jif3D::platform::drand48);
        std::sort(Frequencies.begin(), Frequencies.end(), std::greater<double>());
        std::generate_n(XCoord.begin(), nstat, jif3D::platform::drand48);
        std::generate_n(YCoord.begin(), nstat, jif3D::platform::drand48);
        std::generate_n(ZCoord.begin(), nstat, jif3D::platform::drand48);
        std::generate_n(Impedances.begin(), ndata, jif3D::platform::drand48);
        std::generate_n(Error.begin(), ndata, jif3D::platform::drand48);
        std::generate_n(C.begin(), nstat * 4, jif3D::platform::drand48);
        for (size_t i = 1; i < ndata; i += 2)
          {
            Error(i) = Error(i - 1);
          }

      }

    BOOST_AUTO_TEST_CASE (read_write_netcdf_test)
      {

        std::vector<double> Frequencies;
        std::vector<double> XCoord, YCoord, ZCoord, C;
        jif3D::rvec Impedances, Error;

        GenerateData(Frequencies, XCoord, YCoord, ZCoord, Impedances, Error, C);
        const size_t nfreq = Frequencies.size();
        const size_t nstat = XCoord.size();
        const size_t ndata = nfreq * nstat * 8;
        const std::string filename("imp.nc");
        jif3D::WriteImpedancesToNetCDF(filename, Frequencies, XCoord, YCoord, ZCoord,
            Impedances, Error, C);

        std::vector<double> ReadFrequencies;
        std::vector<double> ReadXCoord, ReadYCoord, ReadZCoord, ReadC;
        jif3D::rvec ReadImpedances, ReadError;
        jif3D::ReadImpedancesFromNetCDF(filename, ReadFrequencies, ReadXCoord, ReadYCoord,
            ReadZCoord, ReadImpedances, ReadError, ReadC);
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
        for (size_t i = 0; i < nstat * 4; ++i)
          {
            BOOST_CHECK_CLOSE(C[i], ReadC[i], 0.001);

          }
        for (size_t i = 0; i < ndata; ++i)
          {
            BOOST_CHECK_CLOSE(Impedances(i), ReadImpedances(i), 0.001);
            BOOST_CHECK_CLOSE(Error(i), ReadError(i), 0.001);
          }
      }

    BOOST_AUTO_TEST_CASE (read_write_tipper_netcdf_test)
      {

        std::vector<double> Frequencies;
        std::vector<double> XCoord, YCoord, ZCoord;
        jif3D::rvec Tipper, Error;

        const size_t nstat = 11;
        const size_t nfreq = 5;
        Frequencies.resize(nfreq);
        XCoord.resize(nstat);
        YCoord.resize(nstat);
        ZCoord.resize(nstat);
        const size_t ndata = nfreq * nstat * 4;
        Tipper.resize(ndata), Error.resize(ndata);

        std::generate_n(Frequencies.begin(), nfreq, jif3D::platform::drand48);
        std::sort(Frequencies.begin(), Frequencies.end(), std::greater<double>());
        std::generate_n(XCoord.begin(), nstat, jif3D::platform::drand48);
        std::generate_n(YCoord.begin(), nstat, jif3D::platform::drand48);
        std::generate_n(ZCoord.begin(), nstat, jif3D::platform::drand48);
        std::generate_n(Tipper.begin(), ndata, jif3D::platform::drand48);
        std::generate_n(Error.begin(), ndata, jif3D::platform::drand48);

        for (size_t i = 1; i < ndata; i += 2)
          {
            Error(i) = Error(i - 1);
          }

        const std::string filename("tipper.nc");
        jif3D::WriteTipperToNetCDF(filename, Frequencies, XCoord, YCoord, ZCoord, Tipper,
            Error);

        std::vector<double> ReadFrequencies;
        std::vector<double> ReadXCoord, ReadYCoord, ReadZCoord;
        jif3D::rvec ReadTipper, ReadError;
        jif3D::ReadTipperFromNetCDF(filename, ReadFrequencies, ReadXCoord, ReadYCoord,
            ReadZCoord, ReadTipper, ReadError);
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

        for (size_t i = 0; i < ndata; ++i)
          {
            BOOST_CHECK_CLOSE(Tipper(i), Tipper(i), 0.001);
            BOOST_CHECK_CLOSE(Error(i), ReadError(i), 0.001);
          }
      }

    BOOST_AUTO_TEST_CASE (read_write_ModEM_test)
      {

        std::vector<double> Frequencies;
        std::vector<double> XCoord, YCoord, ZCoord, C;
        jif3D::rvec Impedances, Error;

        GenerateData(Frequencies, XCoord, YCoord, ZCoord, Impedances, Error, C);
        const size_t nfreq = Frequencies.size();
        const size_t nstat = XCoord.size();
        const size_t ndata = nfreq * nstat * 8;
        const std::string filename("imp.modem");
        jif3D::WriteImpedancesToModEM(filename, Frequencies, XCoord, YCoord, ZCoord,
            Impedances, Error);

        std::vector<double> ReadFrequencies;
        std::vector<double> ReadXCoord, ReadYCoord, ReadZCoord;
        jif3D::rvec ReadImpedances, ReadError;
        jif3D::ReadImpedancesFromModEM(filename, ReadFrequencies, ReadXCoord, ReadYCoord,
            ReadZCoord, ReadImpedances, ReadError);
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

        for (size_t i = 0; i < ndata; ++i)
          {
            BOOST_CHECK_CLOSE(Impedances(i), ReadImpedances(i), 0.001);
            BOOST_CHECK_CLOSE(Error(i), ReadError(i), 0.001);
          }
      }

    BOOST_AUTO_TEST_CASE (read_write_mtt_test)
      {

        std::vector<double> Frequencies;
        std::vector<double> XCoord, YCoord, ZCoord, C;
        jif3D::rvec Impedances, Error, Tipper, TipErr;

        GenerateData(Frequencies, XCoord, YCoord, ZCoord, Impedances, Error, C, 1);
        Tipper.resize(Impedances.size()/2);
        TipErr.resize(Impedances.size()/2);
        std::generate(Tipper.begin(), Tipper.end(), jif3D::platform::drand48);
        std::generate(TipErr.begin(), TipErr.end(), jif3D::platform::drand48);

        const size_t nfreq = Frequencies.size();
        const size_t nstat = XCoord.size();
        const size_t ndata = nfreq * nstat * 8;
        const std::string filename("imp.mtt");
        jif3D::WriteImpedancesToMtt(filename, Frequencies, Impedances, Error, Tipper,
            TipErr);

        std::vector<double> ReadFrequencies;
        jif3D::rvec ReadImpedances, ReadError, ReadTipper, ReadTipErr;
        jif3D::ReadImpedancesFromMTT(filename + "0.mtt", ReadFrequencies, ReadImpedances,
            ReadError, ReadTipper, ReadTipErr);
        for (size_t i = 0; i < nfreq; ++i)
          {
            BOOST_CHECK_CLOSE(Frequencies[i], ReadFrequencies[i], 0.001);
          }

        for (size_t i = 0; i < ndata; ++i)
          {
            BOOST_CHECK_CLOSE(Impedances(i), ReadImpedances(i), 0.001);
            BOOST_CHECK_CLOSE(Error(i), ReadError(i), 0.001);
          }
        for (size_t i = 0; i < ReadTipper.size(); ++i)
          {
            BOOST_CHECK_CLOSE(Tipper(i), Tipper(i), 0.001);
            BOOST_CHECK_CLOSE(TipErr(i), TipErr(i), 0.001);
          }
      }

    BOOST_AUTO_TEST_CASE (read_write_J_test)
      {

        std::vector<double> MttFrequencies;
        jif3D::rvec MttImpedances, MttError, MttTipper, MttTipErr;

        std::vector<double> JFrequencies;
        jif3D::rvec JImpedances, JError;
        double XC, YC, ZC;
        jif3D::ReadImpedancesFromJ("testJ.j", JFrequencies, XC, YC, ZC, JImpedances,
            JError);
        jif3D::ReadImpedancesFromMTT("testJ.mtt", MttFrequencies, MttImpedances,
            MttError,  MttTipper, MttTipErr);
        for (size_t i = 0; i < MttFrequencies.size(); ++i)
          {
            BOOST_CHECK_CLOSE(MttFrequencies[i], JFrequencies[i], 0.001);
          }

        for (size_t i = 0; i < MttImpedances.size(); ++i)
          {
            BOOST_CHECK_CLOSE(MttImpedances(i), JImpedances(i), 0.001);
            BOOST_CHECK_CLOSE(MttError(i), JError(i), 0.001);
          }
      }
    BOOST_AUTO_TEST_SUITE_END()
