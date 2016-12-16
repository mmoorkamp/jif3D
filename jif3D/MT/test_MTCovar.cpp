//============================================================================
// Name        : test_ReadWriteX3D.cpp
// Author      : Feb 17, 2009
// Version     :
// Copyright   : 2009, mmoorkamp
//============================================================================

#define BOOST_TEST_MODULE OneDMTCalculator test
#define BOOST_TEST_MAIN ...
#include "../Global/Jif3DTesting.h"
#include <boost/test/included/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

#include "MTCovar.h"
#include "../Global/FatalException.h"
#include "../Global/Jif3DPlatformHelper.h"

BOOST_AUTO_TEST_SUITE( MTCovar_Suite )

    BOOST_AUTO_TEST_CASE (size_test)
      {
        jif3D::covmat TooBig(5, 5);
        jif3D::covmat TooSmall(2, 2);
        jif3D::covmat RightSize(4, 4);
        BOOST_CHECK_THROW(jif3D::RotateMTCovar(1.0, TooBig), jif3D::FatalException);
        BOOST_CHECK_THROW(jif3D::RotateMTCovar(1.0, TooSmall), jif3D::FatalException);
        BOOST_CHECK_NO_THROW(jif3D::RotateMTCovar(1.0, RightSize));
      }

    BOOST_AUTO_TEST_CASE (identity_test)
      {
        jif3D::covmat RandMat(4, 4);
        jif3D::platform::srand48((unsigned int) time(nullptr));
        for (size_t i = 0; i < 4; ++i)
          for (size_t j = i; j < 4; ++j)
            RandMat(i, j) = jif3D::platform::drand48();

        jif3D::covmat CompMat = jif3D::RotateMTCovar(0.0, RandMat);
        for (size_t i = 0; i < 4; ++i)
          for (size_t j = 0; j < 4; ++j)
            {
              BOOST_CHECK_CLOSE(RandMat(i, j), CompMat(i, j), 0.01);
            }

      }
    BOOST_AUTO_TEST_SUITE_END()
