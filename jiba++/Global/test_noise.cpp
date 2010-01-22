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
#include <boost/bind.hpp>

using namespace boost::accumulators;

BOOST_AUTO_TEST_SUITE( Noise_Test_Suite )

BOOST_AUTO_TEST_CASE  (addnoise)
    {
      const size_t ndata = 1e6;
      jiba::rvec Data(ndata);
      std::fill_n(Data.begin(),ndata,1.0);
      jiba::rvec Noisy(Data);
      //test relative noise level
      const double relerror = 0.1;
      jiba::AddNoise(Noisy,relerror,0.0);
      jiba::rvec Diff(Data-Noisy);

      accumulator_set<double, stats<tag::mean,tag::variance > > acc;
      std::for_each( Diff.begin(), Diff.end(), boost::bind<void>( boost::ref(acc), _1 ) );
      BOOST_CHECK(mean(acc) < 0.001);
      BOOST_CHECK_CLOSE(variance(acc),relerror*relerror,0.2);

      //now check absolute noise level
      Noisy = Data;
      const double abserror = 1;
      jiba::AddNoise(Noisy,0.0,abserror);
      Diff = Data - Noisy;
      accumulator_set<double, stats<tag::mean,tag::variance > > acc2;
      std::for_each( Diff.begin(), Diff.end(), boost::bind<void>( boost::ref(acc2), _1 ) );
      BOOST_CHECK(mean(acc2) < 0.001);
      BOOST_CHECK_CLOSE(variance(acc2),abserror*abserror,0.2);
    }
  BOOST_AUTO_TEST_SUITE_END()
