//============================================================================
// Name        : test_interpolate.cpp
// Author      : Max Moorkamp
// Version     :
// Copyright   : 2008, MM
//============================================================================

//check the quadratic interpolation routine
#define BOOST_TEST_MODULE Interpolation test
#define BOOST_TEST_MAIN ...
#include <boost/test/included/unit_test.hpp>
#include <stdlib.h>
#include "Interpolate.h"
#include <boost/test/floating_point_comparison.hpp>

BOOST_AUTO_TEST_SUITE( Interpolation_Test_Suite )

BOOST_AUTO_TEST_CASE(basic_quad_inter)
      {
       //we try a number of different random quadratic functions
        const size_t tries = 20;
        for (size_t i = 0; i < tries; ++i)
          {
            //generate some coefficients
            double a = drand48();
            double b = drand48();
            double c = drand48();
            //calculate the true minimum
            double min = -b / (2.0 * a);
            //evaluate the function at a few points
            double x1 = min - double(rand());
            double y1 = a * x1 * x1 + b * x1 + c;
            double x2 = min + double(rand());
            double y2 = a * x2 * x2 + b * x2 + c;
            double x3 = min + double(rand());
            double y3 = a * x3 * x3 + b * x3 + c;
            //interpolate the minimum
            double zero = jiba::QuadraticInterpolation(x1, y1, x2, y2, x3, y3);
            //check that it matches the true minimum
            BOOST_CHECK_CLOSE(zero,min,0.01);
          }
      }
    BOOST_AUTO_TEST_SUITE_END()
