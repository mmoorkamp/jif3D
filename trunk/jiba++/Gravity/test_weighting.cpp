//============================================================================
// Name        : test_weighting.cpp
// Author      : Sep 18, 2008
// Version     :
// Copyright   : 2008, mmoorkamp
//============================================================================


#define BOOST_TEST_MODULE GravityDepthWeighting test
#define BOOST_TEST_MAIN ...
#include <boost/test/included/unit_test.hpp>
#include <stdlib.h>
#include <time.h>
#include "DepthWeighting.h"
#include <boost/test/floating_point_comparison.hpp>

BOOST_AUTO_TEST_SUITE( Gravity_DepthWeighting_Suite )

BOOST_AUTO_TEST_CASE    (weightingterm_test)
      {
        //this value is arbitrary
        const double z0 = acos(-1.0);
        const double exponent = -2.0;
        const double z = 0.4;
        double value1 = jiba::WeightingTerm(exponent)(z,z0);
        const double delta = std::numeric_limits<float>::epsilon();
        double value2 = jiba::WeightingTerm(exponent)(z+2.0*delta,z0);
        double deriv = jiba::WeightingTerm(exponent).deriv(z+delta,z0);
        BOOST_CHECK_CLOSE(deriv, (value2-value1)/(2.0*delta), 0.1);
      }

    BOOST_AUTO_TEST_CASE(fitz0_test)
      {
        jiba::ThreeDGravityModel Model;
        const size_t nz =20;
        //we only need the model to make FitZ0 happy
        //all it does is determine the starting index
        //we are only interested in z dependence
        //so we have only 1 cell in x and y
        Model.SetXCellSizes().resize(boost::extents[1]);
        Model.SetYCellSizes().resize(boost::extents[1]);
        Model.SetZCellSizes().resize(boost::extents[nz]);
        Model.SetXCellSizes()[0] = 10.0;
        Model.SetYCellSizes()[0] = 10.0;
        //we make up some "sensitivities"
        jiba::rvec PseudoSens(nz);
        //we want to get back this number
        const double z0 = 237.0;
        //this is the exponent we use to fit the gravity data
        const double exponent = -3.0;
        double currdepth = 0.0;
        //fill cell sizes and "sensitivities" with well defined values
        for (size_t i = 0; i < nz; ++i)
          {
            Model.SetZCellSizes()[i] = nz * i;
            currdepth += Model.GetZCellSizes()[i];
            PseudoSens[i] = jiba::WeightingTerm(exponent)(currdepth,z0);
          }
        //FitZ0 should find z0 specified above
        double foundz0 = jiba::FitZ0(PseudoSens,Model,jiba::WeightingTerm(exponent));
        //in reality the gradient is small at the end, so we can have 1% difference
        BOOST_CHECK_CLOSE(z0, foundz0, 1);
      }

    BOOST_AUTO_TEST_SUITE_END()
