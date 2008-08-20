//============================================================================
// Name        : test_interpolate.cpp
// Author      : Max Moorkamp
// Version     :
// Copyright   : 2008, MM
//============================================================================

#define BOOST_TEST_MODULE Interpolation test
#define BOOST_TEST_MAIN ...
#include <boost/test/included/unit_test.hpp>
#include <stdlib.h>
#include "Interpolate.h"
#include <boost/test/floating_point_comparison.hpp>

BOOST_AUTO_TEST_SUITE( Interpolation_Test_Suite )

BOOST_AUTO_TEST_CASE(basic_quad_inter)
  {
    double a = double(rand());
    double b = double(rand());
    double c = double(rand());
    double min = -b/(2.0 * a);
    double x1 = min - 10.0;
    double y1 = a*x1*x1 + b * x1 +c;
    double x2 = min -5;
    double y2 = a*x2*x2 + b * x2 +c;
    double x3 = min + 10;
    double y3 = a*x3*x3 + b * x3 +c;
    double zero = jiba::QuadraticInterpolation(x1,y1,x2,y2,x3,y3);
    BOOST_CHECK_CLOSE(zero,min,std::numeric_limits<float>::epsilon());
  }
BOOST_AUTO_TEST_SUITE_END()
