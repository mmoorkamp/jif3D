//============================================================================
// Name        : test_interpolate.cpp
// Author      : Max Moorkamp
// Version     :
// Copyright   : 2008, MM
//============================================================================

//check the quadratic interpolation routine
#define BOOST_TEST_MODULE Interpolation test
#define BOOST_TEST_MAIN ...

#include <stdlib.h>
#include "Noise.h"
#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/included/unit_test.hpp>
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/accumulators/statistics/variance.hpp>

using namespace boost::accumulators;

BOOST_AUTO_TEST_SUITE( Noise_Test_Suite )

    BOOST_AUTO_TEST_CASE(relnoise)
      {
        const size_t ndata = 1e4;
        jif3D::rvec Data(ndata);
        std::fill_n(Data.begin(), ndata, 1.0);
        const double relerror = 0.1;
        jif3D::AddNoise(Data, relerror, 0.0);
        for (double value : Data)
          {
            BOOST_CHECK(value != 1.0);
          }
      }

    BOOST_AUTO_TEST_CASE(absnoise)
      {
        const size_t ndata = 1e4;
        jif3D::rvec Data(ndata);
        std::fill_n(Data.begin(), ndata, 1.0);
        const double abserror = 0.1;
        jif3D::AddNoise(Data, 0.0, abserror);
        for (double value : Data)
          {
            BOOST_CHECK(value != 1.0);
          }
      }

    BOOST_AUTO_TEST_CASE(zero)
      {
        const size_t ndata = 1e4;
        jif3D::rvec Data(ndata);
        std::fill_n(Data.begin(), ndata, 1.0);
        jif3D::AddNoise(Data, 0.0, 0.0);
        for (double value : Data)
          {
            BOOST_CHECK(value == 1.0);
          }
      }

//check that the function to add noise to synthetic data
//produces a Gaussian distribution with the required parameters
//this test can fail occassionally due to the finite sample size
    BOOST_AUTO_TEST_CASE (var)
      {
        //we use a million data points
        const size_t ndata = 1e6;
        jif3D::rvec Data(ndata);
        std::fill_n(Data.begin(), ndata, 1.0);
        jif3D::rvec Noisy(Data);
        //test relative noise level
        const double relerror = 0.1;
        jif3D::AddNoise(Noisy, relerror, 0.0);
        jif3D::rvec Diff(Data - Noisy);

        accumulator_set<double, stats<tag::mean, tag::variance> > acc;
        std::for_each(Diff.begin(), Diff.end(), [&acc] (double d) {acc(d);});
        BOOST_CHECK(mean(acc) < 0.001);
        BOOST_CHECK_CLOSE(variance(acc), relerror * relerror, 0.2);

        //now check absolute noise level
        Noisy = Data;
        const double abserror = 1;
        jif3D::AddNoise(Noisy, 0.0, abserror);
        Diff = Data - Noisy;
        accumulator_set<double, stats<tag::mean, tag::variance> > acc2;
        std::for_each(Diff.begin(), Diff.end(), [&acc2] (double d) {acc2(d);});
        BOOST_CHECK(mean(acc2) < 0.001);
        BOOST_CHECK_CLOSE(variance(acc2), abserror * abserror, 0.2);
      }

    BOOST_AUTO_TEST_SUITE_END()
