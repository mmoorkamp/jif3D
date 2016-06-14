/*
 * test_transform.cpp
 *
 *  Created on: 21 Apr 2016
 *      Author: mm489
 */

#define BOOST_TEST_MODULE MTTransforms test
#define BOOST_TEST_MAIN ...
#include <boost/test/included/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>
//#include "test_common.h"
#include "MTTransforms.h"
#include "MTEquations.h"
#include "../Global/NumUtil.h"

BOOST_AUTO_TEST_SUITE( MT_Transforms_Suite )

    BOOST_AUTO_TEST_CASE(basic_trans_test)
      {
        const double freq = 10.0;
        const double cond = 0.05;
        std::complex<double> Imp1D = jif3D::ImpedanceHalfspace(freq, cond);
        jif3D::ComplexLogTransform Trans;
        jif3D::rvec values(8);
        for (size_t i = 0; i < values.size(); i += 2)
          {
            values(i) = Imp1D.real();
            values(i + 1) = Imp1D.imag();
          }
        jif3D::rvec TValues = jif3D::ApplyTransform(values, Trans);
        BOOST_CHECK_EQUAL(TValues.size(), values.size());
        for (size_t i = 0; i < values.size(); i += 2)
          {
            BOOST_CHECK_CLOSE(1.0 / (jif3D::twopimu * freq) * exp(2.0 * TValues(i)),
                1.0 / cond, 0.01);
            BOOST_CHECK_CLOSE(
                180.0 / boost::math::constants::pi<double>() * TValues(i + 1), 45.0,
                0.01);
          }
      }

    BOOST_AUTO_TEST_CASE(finitediff_test)
      {
        jif3D::ComplexLogTransform CLTrans;
        jif3D::rvec TestInput(2);
        std::generate(TestInput.begin(), TestInput.end(), rand);

        jif3D::rvec CL(CLTrans.Transform(TestInput));
        jif3D::rmat Deriv(CLTrans.Derivative(TestInput));
        const double delta = 0.001;
        for (size_t i = 0; i < TestInput.size(); ++i)
          {
            for (size_t j = 0; j < TestInput.size(); ++j)
              {
                jif3D::rvec DeltaInput(TestInput);
                DeltaInput(i) *= 1.0 + delta;
                double FD = (CLTrans.Transform(DeltaInput)(j) - CL(j))
                    / (delta * TestInput(i));
                BOOST_CHECK_CLOSE(Deriv(j,i), FD, 0.1);
              }
          }
      }
    BOOST_AUTO_TEST_SUITE_END()

