//============================================================================
// Name        : test_wavelet.cpp
// Author      : Max Moorkamp
// Version     :
// Copyright   : 2008, MM
//============================================================================

//test the wavelet forward and inverse transforms
#define BOOST_TEST_MODULE FFT test
#define BOOST_TEST_MAIN ...

#include "VecMat.h"
#include "kissfft.hh"
#include "Jif3DPlatformHelper.h"
#include "Jif3DTesting.h"
//#include <boost/test/tools/floating_point_comparison.hpp>



BOOST_AUTO_TEST_SUITE( FFT_Test_Suite )

    BOOST_AUTO_TEST_CASE (FFT_transform_pair)
      {
        //we check that a transform followed by the inverse
        //gives back the original result for a simple vector

        //construct the vector and fill it with random numbers
        //we use integers, because drand48 can generate very small values
        //that reach the numerical precision and give false alarms
        const size_t length = 64;
        jif3D::cvec Vector(length), Result(length);
        std::generate_n(Vector.begin(), length, drand48);
        //store the original
        jif3D::cvec Original(Vector);
        //do the forward
        kissfft<double> fft(length,false);
        fft.transform(&Vector.data()[0],&Result.data()[0]);
        //and the inverse
        kissfft<double> ifft(length,true);
        ifft.transform(&Result.data()[0],&Vector.data()[0]);
        //check that they are still the same
        for (size_t i = 0; i < length; ++i)
          {
            BOOST_CHECK_CLOSE(Original(i).real(), 1.0/(length) * Vector(i).real(), 1e-3);
            BOOST_CHECK( std::abs(Vector(i).imag()) < 1e-10);
          }
      }

   BOOST_AUTO_TEST_SUITE_END()
