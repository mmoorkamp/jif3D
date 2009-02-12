//============================================================================
// Name        : test_interpolate.cpp
// Author      : Max Moorkamp
// Version     :
// Copyright   : 2008, MM
//============================================================================

#define BOOST_TEST_MODULE NumUtil test
#define BOOST_TEST_MAIN ...
#include <boost/test/included/unit_test.hpp>
#include <stdlib.h>
#include "NumUtil.h"
#include <boost/test/floating_point_comparison.hpp>

BOOST_AUTO_TEST_SUITE( NumUtil_Test_Suite )

BOOST_AUTO_TEST_CASE(test_sign)
    {
      double posdoub = 10.0;
      double negdoub = -10.0;
      int posint = 5;
      BOOST_CHECK( jiba::sign(posdoub) ==1);
      BOOST_CHECK( jiba::sign(negdoub) ==-1);
      BOOST_CHECK( jiba::sign(posint) ==1);
    }

  BOOST_AUTO_TEST_CASE(test_intseq)
    {
      const size_t ntries = 30;
      const int start = 7;
      std::vector<int> Sequence;
      std::generate_n(back_inserter(Sequence), ntries, jiba::IntSequence(start));
      for (size_t i = 0; i < ntries; ++i)
        {
          BOOST_CHECK(Sequence.at(i)== i+start);
        }

    }

  BOOST_AUTO_TEST_CASE(test_absless)
    {
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

BOOST_AUTO_TEST_SUITE_END()
