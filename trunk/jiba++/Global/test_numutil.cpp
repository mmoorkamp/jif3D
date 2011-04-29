//============================================================================
// Name        : test_numutil.cpp
// Author      : Max Moorkamp
// Version     :
// Copyright   : 2008, MM
//============================================================================


//test the small numerical utility functions
#define BOOST_TEST_MODULE NumUtil test
#define BOOST_TEST_MAIN ...
#include <boost/test/included/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <stdlib.h>
#include "NumUtil.h"
#include "NormProd.h"
#include "VecMat.h"

BOOST_AUTO_TEST_SUITE( NumUtil_Test_Suite )

BOOST_AUTO_TEST_CASE  (test_sign)
    {
      //an extremely simple test for the sign function
      double posdoub = 10.0;
      double negdoub = -10.0;
      int posint = 5;
      BOOST_CHECK( jiba::sign(posdoub) ==1);
      BOOST_CHECK( jiba::sign(negdoub) ==-1);
      BOOST_CHECK( jiba::sign(posint) ==1);
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
      std::generate_n(back_inserter(Sequence), ntries, jiba::IntSequence(start));
      for (int i = 0; i < ntries; ++i)
        {
          BOOST_CHECK(Sequence.at(i)== i+start);
        }

    }

  BOOST_AUTO_TEST_CASE(test_absless)
    {
      //test different instances of the absolute less function
      double posdoub = 1.0;
      double negdoub = -10.0;
      BOOST_CHECK( (jiba::absLess<double,double>()(posdoub,negdoub) == true) );
      BOOST_CHECK( ( jiba::absLess<double,double>()(negdoub,posdoub) == false) );
      negdoub = -1.0;
      BOOST_CHECK( (jiba::absLess<double,double>()(posdoub,negdoub) == false) );
      negdoub = -0.5;
      BOOST_CHECK( ( jiba::absLess<double,double>()(posdoub,negdoub) == false) );
      BOOST_CHECK( ( jiba::absLess<double,double>()(negdoub,posdoub) == true) );
    }

  BOOST_AUTO_TEST_CASE(test_normprod)
    {
      //check that our normed scalar product function works correctly
      const size_t nelements = 1010;
      jiba::rvec a(nelements), b(nelements), NormDiag(nelements);
      std::generate_n(a.begin(),nelements, drand48);
      std::generate_n(b.begin(),nelements, drand48);
      std::generate_n(NormDiag.begin(),nelements, drand48);
      double myresult = jiba::NormProd(a, b, NormDiag);
      double uresult = ublas::inner_prod(a, ublas::element_div(b, NormDiag));
      BOOST_CHECK_CLOSE(myresult,uresult,std::numeric_limits<float>::epsilon());
    }

  BOOST_AUTO_TEST_CASE(test_poweroftwo)
    {
      //generate an exponent between 2 and 22
      size_t n = (rand() % 20) + 2;
      size_t ispower1 = std::pow(2,n);
      BOOST_CHECK(jiba::IsPowerOfTwo(ispower1));
      size_t m = n +3;
      size_t ispower2 = std::pow(2,m);
      BOOST_CHECK(jiba::IsPowerOfTwo(ispower2));
      //the sum of the two numbers with different exponents cannot be a power of two
      BOOST_CHECK(!jiba::IsPowerOfTwo(ispower1+ispower2));
    }
  BOOST_AUTO_TEST_SUITE_END()
