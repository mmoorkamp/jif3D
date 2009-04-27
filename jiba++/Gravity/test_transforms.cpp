//============================================================================
// Name        : test_transforms.cpp
// Author      : Apr 22, 2009
// Version     :
// Copyright   : 2009, mmoorkamp
//============================================================================


#define BOOST_TEST_MODULE GravityTransforms test
#define BOOST_TEST_MAIN ...
#include <boost/test/included/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <stdlib.h>
#include <time.h>
#include "GravityTransforms.h"
#include "../Global/NumUtil.h"

BOOST_AUTO_TEST_SUITE( Gravity_Transforms_Suite )

BOOST_AUTO_TEST_CASE(ftginvar_test)
    {
     jiba::FTGInvariant InvarTrans;
     jiba::rvec TestInput(9);
     TestInput(0) = 1.0;
     TestInput(1) = 2.0;
     TestInput(2) = 3.0;
     TestInput(3) = 2.0;
     TestInput(4) = 4.0;
     TestInput(5) = 5.0;
     TestInput(6) = 3.0;
     TestInput(7) = 5.0;
     TestInput(8) = 6.0;
     jiba::rvec Invar(InvarTrans.Transform(TestInput));
     jiba::rmat Deriv(InvarTrans.Derivative(TestInput));
     BOOST_CHECK(Invar.size() == 1);
     BOOST_CHECK(Deriv.size1() == 1);
     BOOST_CHECK(Deriv.size2() == 9);
     BOOST_CHECK_CLOSE(Invar(0),-4.0,std::numeric_limits<float>::epsilon());
     BOOST_CHECK_CLOSE(Deriv(0,0),10.0,std::numeric_limits<float>::epsilon());
     BOOST_CHECK_CLOSE(Deriv(0,1),-4.0,std::numeric_limits<float>::epsilon());
     BOOST_CHECK_CLOSE(Deriv(0,2),-6.0,std::numeric_limits<float>::epsilon());
     BOOST_CHECK_CLOSE(Deriv(0,3),-4.0,std::numeric_limits<float>::epsilon());
     BOOST_CHECK_CLOSE(Deriv(0,4),7.0,std::numeric_limits<float>::epsilon());
     BOOST_CHECK_CLOSE(Deriv(0,5),-10.0,std::numeric_limits<float>::epsilon());
     BOOST_CHECK_CLOSE(Deriv(0,6),-6.0,std::numeric_limits<float>::epsilon());
     BOOST_CHECK_CLOSE(Deriv(0,7),-10.0,std::numeric_limits<float>::epsilon());
     BOOST_CHECK_CLOSE(Deriv(0,8),5.0,std::numeric_limits<float>::epsilon());

    }

BOOST_AUTO_TEST_SUITE_END()
