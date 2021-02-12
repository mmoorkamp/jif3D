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
#include <boost/test/tools/floating_point_comparison.hpp>
#include <random>
#include <string>

#include "../Global/Jif3DPlatformHelper.h"
#include "ReadWriteImpedances.h"
#include "X3DModel.h"
BOOST_AUTO_TEST_SUITE (ReadWriteImpedances_Suite)

    std::string random_string(std::string::size_type length)
      {
        static auto &chrs = "0123456789"
            "abcdefghijklmnopqrstuvwxyz"
            "ABCDEFGHIJKLMNOPQRSTUVWXYZ";

        thread_local static std::mt19937 rg
          { std::random_device
            { }() };
        thread_local static std::uniform_int_distribution<std::string::size_type> pick(0,
            sizeof(chrs) - 2);

        std::string s;

        s.reserve(length);

        while (length--)
          s += chrs[pick(rg)];

        return s;
      }

    void GenerateData(std::vector<double> &Frequencies, std::vector<double> &XCoord,
        std::vector<double> &YCoord, std::vector<double> &ZCoord,
        std::vector<double> &Impedances, std::vector<double> &Error,
        std::vector<double> &C, size_t nstat = 7)
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
            Error.at(i) = Error.at(i - 1);
          }

      }

    BOOST_AUTO_TEST_CASE (read_write_netcdf_test)
      {

        std::vector<double> Frequencies;
        std::vector<double> XCoord, YCoord, ZCoord, C;
        std::vector<double> Impedances, Error;

        GenerateData(Frequencies, XCoord, YCoord, ZCoord, Impedances, Error, C);
        const size_t nfreq = Frequencies.size();
        const size_t nstat = XCoord.size();
        const size_t ndata = nfreq * nstat * 8;
        std::vector<std::string> Names(nstat);
        std::generate_n(Names.begin(), nstat, []()
          { return random_string(19);});
        const std::string filename("imp.nc");
        jif3D::WriteImpedancesToNetCDF(filename, Frequencies, XCoord, YCoord, ZCoord,
            Impedances, Error, C, Names);

        std::vector<double> ReadFrequencies;
        std::vector<double> ReadXCoord, ReadYCoord, ReadZCoord, ReadC;
        std::vector<double> ReadImpedances, ReadError;
        std::vector<std::string> ReadNames;
        jif3D::ReadImpedancesFromNetCDF(filename, ReadFrequencies, ReadXCoord, ReadYCoord,
            ReadZCoord, ReadImpedances, ReadError, ReadC, ReadNames);
        for (size_t i = 0; i < nfreq; ++i)
          {
            BOOST_CHECK_CLOSE(Frequencies[i], ReadFrequencies[i], 0.001);
          }
        for (size_t i = 0; i < nstat; ++i)
          {
            BOOST_CHECK_CLOSE(XCoord[i], ReadXCoord[i], 0.001);
            BOOST_CHECK_CLOSE(YCoord[i], ReadYCoord[i], 0.001);
            BOOST_CHECK_CLOSE(ZCoord[i], ReadZCoord[i], 0.001);
            BOOST_TEST(Names.at(i) == ReadNames.at(i));
          }
        for (size_t i = 0; i < nstat * 4; ++i)
          {
            BOOST_CHECK_CLOSE(C[i], ReadC[i], 0.001);

          }
        for (size_t i = 0; i < ndata; ++i)
          {
            BOOST_CHECK_CLOSE(Impedances.at(i), ReadImpedances.at(i), 0.001);
            BOOST_CHECK_CLOSE(Error.at(i), ReadError.at(i), 0.001);
          }
      }

    BOOST_AUTO_TEST_CASE (read_write_tipper_netcdf_test)
      {

        const size_t nstat = 11;
        const size_t nfreq = 5;
        const size_t ndata = nfreq * nstat * 4;

        std::vector<double> Frequencies(nfreq);
        std::vector<double> XCoord(nstat), YCoord(nstat), ZCoord(nstat);
        std::vector<double> Tipper(ndata), Error(ndata);
        std::vector<int> HxIndex(nstat * nfreq), HyIndex(nstat * nfreq), HzIndex(
            nstat * nfreq);
        std::vector<std::string> Names(nstat);
        std::generate_n(Names.begin(), nstat, []()
          { return random_string(19);});

        std::generate_n(Frequencies.begin(), nfreq, jif3D::platform::drand48);
        std::sort(Frequencies.begin(), Frequencies.end(), std::greater<double>());
        std::generate_n(XCoord.begin(), nstat, jif3D::platform::drand48);
        std::generate_n(YCoord.begin(), nstat, jif3D::platform::drand48);
        std::generate_n(ZCoord.begin(), nstat, jif3D::platform::drand48);
        std::generate_n(Tipper.begin(), ndata, jif3D::platform::drand48);
        std::generate_n(Error.begin(), ndata, jif3D::platform::drand48);
        std::generate_n(HxIndex.begin(), nstat * nfreq, std::rand);
        std::generate_n(HyIndex.begin(), nstat * nfreq, std::rand);
        std::generate_n(HzIndex.begin(), nstat * nfreq, std::rand);

        for (size_t i = 1; i < ndata; i += 2)
          {
            Error.at(i) = Error.at(i - 1);
          }

        const std::string filename("tipper.nc");
        jif3D::WriteTipperToNetCDF(filename, Frequencies, XCoord, YCoord, ZCoord, HxIndex,
            HyIndex, HzIndex, Tipper, Error, Names);

        std::vector<double> ReadFrequencies;
        std::vector<double> ReadXCoord, ReadYCoord, ReadZCoord;
        std::vector<double> ReadTipper, ReadError;
        std::vector<int> ReadHxIndex, ReadHyIndex, ReadHzIndex;
        std::vector<std::string> ReadNames;

        jif3D::ReadTipperFromNetCDF(filename, ReadFrequencies, ReadXCoord, ReadYCoord,
            ReadZCoord, ReadHxIndex, ReadHyIndex, ReadHzIndex, ReadTipper, ReadError,
            ReadNames);
        for (size_t i = 0; i < nfreq; ++i)
          {
            BOOST_CHECK_CLOSE(Frequencies[i], ReadFrequencies[i], 0.001);
          }
        for (size_t i = 0; i < nstat; ++i)
          {
            BOOST_CHECK_CLOSE(XCoord[i], ReadXCoord[i], 0.001);
            BOOST_CHECK_CLOSE(YCoord[i], ReadYCoord[i], 0.001);
            BOOST_CHECK_CLOSE(ZCoord[i], ReadZCoord[i], 0.001);
            BOOST_CHECK_EQUAL(HxIndex[i], ReadHxIndex[i]);
            BOOST_CHECK_EQUAL(HyIndex[i], ReadHyIndex[i]);
            BOOST_CHECK_EQUAL(HzIndex[i], ReadHzIndex[i]);
            BOOST_TEST(Names.at(i) == ReadNames.at(i));

          }

        for (size_t i = 0; i < ndata; ++i)
          {
            BOOST_CHECK_CLOSE(Tipper.at(i), Tipper.at(i), 0.001);
            BOOST_CHECK_CLOSE(Error.at(i), ReadError.at(i), 0.001);
          }
      }

    BOOST_AUTO_TEST_CASE (read_write_ModEM_test)
      {

        std::vector<double> Frequencies;
        std::vector<double> XCoord, YCoord, ZCoord, C;
        std::vector<double> Impedances, Error;

        GenerateData(Frequencies, XCoord, YCoord, ZCoord, Impedances, Error, C);
        const size_t nfreq = Frequencies.size();
        const size_t nstat = XCoord.size();
        const size_t ndata = nfreq * nstat * 8;
        const std::string filename("imp.modem");
        std::vector<std::string> Names(nstat);
        std::generate(Names.begin(), Names.end(), []()
          { return random_string(8);});
        jif3D::WriteImpedancesToModEM(filename, Frequencies, XCoord, YCoord, ZCoord,
            Impedances, Error, Names);

        std::vector<double> ReadFrequencies;
        std::vector<double> ReadXCoord, ReadYCoord, ReadZCoord;
        std::vector<double> ReadImpedances, ReadError;
        std::vector < std::string > ReadNames;
        jif3D::ReadImpedancesFromModEM(filename, ReadFrequencies, ReadXCoord, ReadYCoord,
            ReadZCoord, ReadImpedances, ReadError, ReadNames);
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
            BOOST_CHECK_CLOSE(Impedances.at(i), ReadImpedances.at(i), 0.001);
            BOOST_CHECK_CLOSE(Error.at(i), ReadError.at(i), 0.001);
          }
      }

    BOOST_AUTO_TEST_CASE (read_write_mtt_test)
      {

        std::vector<double> Frequencies;
        std::vector<double> XCoord, YCoord, ZCoord, C;
        std::vector<double> Impedances, Error, Tipper, TipErr;

        GenerateData(Frequencies, XCoord, YCoord, ZCoord, Impedances, Error, C, 1);
        Tipper.resize(Impedances.size() / 2);
        TipErr.resize(Impedances.size() / 2);
        std::generate(Tipper.begin(), Tipper.end(), jif3D::platform::drand48);
        std::generate(TipErr.begin(), TipErr.end(), jif3D::platform::drand48);

        const size_t nfreq = Frequencies.size();
        const size_t nstat = XCoord.size();
        const size_t ndata = nfreq * nstat * 8;
        const std::string filename("imp.mtt");
        jif3D::WriteImpedancesToMtt(filename, Frequencies, Impedances, Error, Tipper,
            TipErr);

        std::vector<double> ReadFrequencies;
        std::vector<double> ReadImpedances, ReadError, ReadTipper, ReadTipErr;
        jif3D::ReadImpedancesFromMTT(filename + "0.mtt", ReadFrequencies, ReadImpedances,
            ReadError, ReadTipper, ReadTipErr);
        for (size_t i = 0; i < nfreq; ++i)
          {
            BOOST_CHECK_CLOSE(Frequencies[i], ReadFrequencies[i], 0.001);
          }

        for (size_t i = 0; i < ndata; ++i)
          {
            BOOST_CHECK_CLOSE(Impedances.at(i), ReadImpedances.at(i), 0.001);
            BOOST_CHECK_CLOSE(Error.at(i), ReadError.at(i), 0.001);
          }
        for (size_t i = 0; i < ReadTipper.size(); ++i)
          {
            BOOST_CHECK_CLOSE(Tipper.at(i), Tipper.at(i), 0.001);
            BOOST_CHECK_CLOSE(TipErr.at(i), TipErr.at(i), 0.001);
          }
      }

    BOOST_AUTO_TEST_CASE (read_write_J_test)
      {

        std::vector<double> MttFrequencies;
        std::vector<double> MttImpedances, MttError, MttTipper, MttTipErr;

        std::vector<double> JFrequencies;
        std::vector<double> JImpedances, JError, JTipper, JTipErr;
        double XC, YC, ZC;
        jif3D::ReadImpedancesFromJ("testJ.j", JFrequencies, XC, YC, ZC, JImpedances,
            JError);
        jif3D::ReadImpedancesFromMTT("testJ.mtt", MttFrequencies, MttImpedances, MttError,
            MttTipper, MttTipErr);
        for (size_t i = 0; i < MttFrequencies.size(); ++i)
          {
            BOOST_CHECK_CLOSE(MttFrequencies[i], JFrequencies[i], 0.001);
          }

        for (size_t i = 0; i < MttImpedances.size(); ++i)
          {
            BOOST_CHECK_CLOSE(MttImpedances.at(i), JImpedances.at(i), 0.001);
            BOOST_CHECK_CLOSE(MttError.at(i), JError.at(i), 0.001);
          }
      }
    BOOST_AUTO_TEST_SUITE_END()
