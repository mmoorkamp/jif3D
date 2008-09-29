//============================================================================
// Name        : test_wavelet.cpp
// Author      : Max Moorkamp
// Version     :
// Copyright   : 2008, MM
//============================================================================

#define BOOST_TEST_MODULE Interpolation test
#define BOOST_TEST_MAIN ...
#include <boost/test/included/unit_test.hpp>
#include <stdlib.h>
#include "Wavelet.h"
#include <boost/test/floating_point_comparison.hpp>

BOOST_AUTO_TEST_SUITE( Interpolation_Test_Suite )

BOOST_AUTO_TEST_CASE    (wavelet_transform_pair)
      {
        const size_t length = 64;
        jiba::rvec Vector(length);
        std::generate_n(Vector.begin(),length,rand);
        Vector(length/2) = 1.0;
        jiba::rvec Original(Vector);
        jiba::WaveletTransform(Vector);
        jiba::InvWaveletTransform(Vector);
        for (size_t i = 0; i < length; ++i)
          {
            BOOST_CHECK_CLOSE(Original(i),Vector(i),1e-3);
          }
      }
    BOOST_AUTO_TEST_SUITE_END()
