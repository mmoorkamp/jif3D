//============================================================================
// Name        : test_inversion.h
// Author      : Max Moorkamp
// Version     :
// Copyright   : 2008, MM
//============================================================================

#define BOOST_TEST_MODULE Distributor test
#define BOOST_TEST_MAIN ...
#include <boost/test/included/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <stdlib.h>
#include "ModelDistributor.h"

BOOST_AUTO_TEST_SUITE( Distributor_Test_Suite )

BOOST_AUTO_TEST_CASE  (basic_copy_test)
    {
      jiba::ModelDistributor Distributor;
      Distributor.AddTransformer(boost::shared_ptr<jiba::VectorTransform>(new jiba::CopyTransform()));
      Distributor.AddTransformer(boost::shared_ptr<jiba::VectorTransform>(new jiba::CopyTransform()));
      jiba::rvec Input(50);
      std::generate(Input.begin(),Input.end(),rand);

      jiba::rvec Output1(Distributor(Input,0));
      jiba::rvec Output2(Distributor(Input,1));

      BOOST_CHECK(std::equal(Input.begin(),Input.end(),Output1.begin()));
      BOOST_CHECK(std::equal(Input.begin(),Input.end(),Output2.begin()));
    }

  BOOST_AUTO_TEST_SUITE_END()
