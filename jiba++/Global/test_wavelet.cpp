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

BOOST_AUTO_TEST_SUITE( Wavelet_Test_Suite )

BOOST_AUTO_TEST_CASE (wavelet_transform_pair)
      {
        //we check that a transform followed by the inverse
        //gives back the original result
        const size_t length = 64;
        jiba::rvec Vector(length);
        std::generate_n(Vector.begin(), length, rand);

        jiba::rvec Original(Vector);
        jiba::WaveletTransform(Vector);
        jiba::InvWaveletTransform(Vector);
        for (size_t i = 0; i < length; ++i)
          {
            BOOST_CHECK_CLOSE(Original(i),Vector(i),1e-3);
          }
      }
    BOOST_AUTO_TEST_CASE (wavelet_output)
      {
        const size_t length = 1024;
        jiba::rvec Vector(length);
        std::fill_n(Vector.begin(),length,0.0);
        Vector(4) = 1.0;
        jiba::InvWaveletTransform(Vector);
        std::ofstream outfile("wavelet.out");
        std::copy(Vector.begin(),Vector.end(),std::ostream_iterator<double>(outfile,"\n"));
      }
    BOOST_AUTO_TEST_SUITE_END()
