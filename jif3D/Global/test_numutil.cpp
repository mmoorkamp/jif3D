//============================================================================
// Name        : test_numutil.cpp
// Author      : Max Moorkamp
// Version     :
// Copyright   : 2008, MM
//============================================================================

//test the small numerical utility functions
#define BOOST_TEST_MODULE NumUtil test
#define BOOST_TEST_MAIN ...

#include "Jif3DTesting.h"
#include <boost/test/floating_point_comparison.hpp>
#include <stdlib.h>
#include "NumUtil.h"
#include "NormProd.h"
#include "VecMat.h"
#include "Jif3DPlatformHelper.h"

BOOST_AUTO_TEST_SUITE( NumUtil_Test_Suite )

    BOOST_AUTO_TEST_CASE (test_sign)
      {
        //an extremely simple test for the sign function
        double posdoub = 10.0;
        double negdoub = -10.0;
        int posint = 5;
        BOOST_CHECK(jif3D::sign(posdoub) == 1);
        BOOST_CHECK(jif3D::sign(negdoub) == -1);
        BOOST_CHECK(jif3D::sign(posint) == 1);
      }

    BOOST_AUTO_TEST_CASE(test_intseq)
      {
        //test the sequence generator
        //the length of the vector
        const int ntries = 30;
        //the value of the first element
        const int start = 7;
        std::vector<int> Sequence;
        //use the generator and compare the result
        std::generate_n(back_inserter(Sequence), ntries, jif3D::IntSequence(start));
        for (int i = 0; i < ntries; ++i)
          {
            BOOST_CHECK(Sequence.at(i) == i + start);
          }

      }

    BOOST_AUTO_TEST_CASE(test_absless)
      {
        //test different instances of the absolute less function
        double posdoub = 1.0;
        double negdoub = -10.0;
        BOOST_CHECK((jif3D::absLess<double, double>()(posdoub, negdoub) == true));
        BOOST_CHECK((jif3D::absLess<double, double>()(negdoub, posdoub) == false));
        negdoub = -1.0;
        BOOST_CHECK((jif3D::absLess<double, double>()(posdoub, negdoub) == false));
        negdoub = -0.5;
        BOOST_CHECK((jif3D::absLess<double, double>()(posdoub, negdoub) == false));
        BOOST_CHECK((jif3D::absLess<double, double>()(negdoub, posdoub) == true));
      }

    BOOST_AUTO_TEST_CASE(test_normprod)
      {
        //check that our normed scalar product function works correctly
        const size_t nelements = 1010;
        jif3D::rvec a(nelements), b(nelements), NormDiag(nelements);
        std::generate_n(a.begin(), nelements, jif3D::platform::drand48);
        std::generate_n(b.begin(), nelements, jif3D::platform::drand48);
        std::generate_n(NormDiag.begin(), nelements, jif3D::platform::drand48);
        double myresult = jif3D::NormProd(a, b, NormDiag);
        double uresult = ublas::inner_prod(a, ublas::element_div(b, NormDiag));
        BOOST_CHECK_CLOSE(myresult, uresult, std::numeric_limits<float>::epsilon());
      }

    BOOST_AUTO_TEST_CASE(test_poweroftwo)
      {
        //generate an exponent between 2 and 22
        size_t n = (rand() % 20) + 2;
        size_t ispower1 = std::pow(2, n);
        BOOST_CHECK(jif3D::IsPowerOfTwo(ispower1));
        size_t m = n + 3;
        size_t ispower2 = std::pow(2, m);
        BOOST_CHECK(jif3D::IsPowerOfTwo(ispower2));
        //the sum of the two numbers with different exponents cannot be a power of two
        BOOST_CHECK(!jif3D::IsPowerOfTwo(ispower1 + ispower2));
      }

    BOOST_AUTO_TEST_CASE(test_roughlyEqual)
      {
        jif3D::platform::srand48(time(NULL));
        const size_t nruns = 10;
        for (size_t i = 0; i < nruns; ++i)
          {
            const double delta = 1e-6 + jif3D::platform::drand48();
            jif3D::roughlyEqual<double,double> Comp(delta);
            double number1 = jif3D::platform::drand48();
            double number2 = number1 + delta / 2.0;
            double number3 = number1 + delta * 2.0;
            BOOST_CHECK(Comp(number1, number2));
            BOOST_CHECK(!Comp(number1, number3));
          }
      }
    BOOST_AUTO_TEST_SUITE_END()
