//============================================================================
// Name        : test_wavelet.cpp
// Author      : Max Moorkamp
// Version     :
// Copyright   : 2008, MM
//============================================================================


//test the wavelet forward and inverse transforms
#define BOOST_TEST_MODULE Wavelet test
#define BOOST_TEST_MAIN ...

#include "Wavelet.h"
#include "Jif3DPlatformHelper.h"
#include <stdlib.h>
#include <boost/multi_array.hpp>
#include "Jif3DTesting.h"
#include <boost/test/floating_point_comparison.hpp>

BOOST_AUTO_TEST_SUITE( Wavelet_Test_Suite )

BOOST_AUTO_TEST_CASE (wavelet_transform_pair)
      {
        //we check that a transform followed by the inverse
        //gives back the original result for a simple vector

        //construct the vector and fill it with random numbers
        //we use integers, because drand48 can generate very small values
        //that reach the numerical precision and give false alarms
        const size_t length = 64;
        jif3D::rvec Vector(length);
        std::generate_n(Vector.begin(), length, rand);
        //store the original
        jif3D::rvec Original(Vector);
        //do the forward
        jif3D::WaveletTransform(Vector);
        //and the inverse
        jif3D::InvWaveletTransform(Vector);
        //check that they are still the same
        for (size_t i = 0; i < length; ++i)
          {
            BOOST_CHECK_CLOSE(Original(i),Vector(i),1e-3);
          }
      }


    BOOST_AUTO_TEST_CASE (multi_dim_1dcomp)
      {
        //we check the forward inverse pair
        //for the multidimensional transform
        //but with an effectively 1D vector
        //and compare it to the "truely" 1D transform

        //first make a 1D vector
        const size_t length = 64;
        jif3D::rvec Vector(length);
        std::generate_n(Vector.begin(), length, jif3D::platform::drand48);
        //then a multi_array that is 1D and copy the vector
        boost::multi_array<double, 1> InArray(boost::extents[length]);
        std::copy(Vector.begin(), Vector.end(), InArray.origin());
        //forward transform both and check that they match
        jif3D::WaveletTransform(Vector);
        jif3D::WaveletTransform(InArray);
        for (size_t i = 0; i < length; ++i)
          {
            BOOST_CHECK_CLOSE(Vector(i),*(InArray.origin()+i),1e-3);
          }
        //inverse transform both and check
        jif3D::InvWaveletTransform(Vector);
        jif3D::InvWaveletTransform(InArray);
        for (size_t i = 0; i < length; ++i)
          {
            BOOST_CHECK_CLOSE(Vector(i),*(InArray.origin()+i),1e-3);
          }
      }

    BOOST_AUTO_TEST_CASE (multi_dim_wavelet_transform_pair)
      {
        //we check that a transform followed by the inverse
        //gives back the original result for 3D arrays
        const size_t length = 8;
        boost::multi_array<double, 3> InArray(
            boost::extents[length][length][length]);
        std::generate_n(InArray.origin(), pow(length, 3), rand);

        boost::multi_array<double, 3> Original(InArray);
        jif3D::WaveletTransform(InArray);
        jif3D::InvWaveletTransform(InArray);
        for (size_t i = 0; i < InArray.num_elements(); ++i)
          {
            BOOST_CHECK_CLOSE(*(Original.origin()+i),*(InArray.origin()+i),1e-3);
          }
      }
    BOOST_AUTO_TEST_SUITE_END()
