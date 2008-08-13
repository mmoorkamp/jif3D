//============================================================================
// Name        : test_inversion.h
// Author      : Max Moorkamp
// Version     :
// Copyright   : 2008, MM
//============================================================================

#define BOOST_TEST_MODULE Inversion test
#define BOOST_TEST_MAIN ...
#include <boost/test/included/unit_test.hpp>
#include <stdlib.h>
#include "MatrixTools.h"
#include <boost/test/floating_point_comparison.hpp>

BOOST_AUTO_TEST_SUITE( Inversion_Test_Suite )

//test that the SVD fulfills the absolute minimum requirement
//calculate the SVD of an Identity matrix
BOOST_AUTO_TEST_CASE(basic_svd_test)
  {
    const size_t msize = 4;
    jiba::rmat TestMatrix(msize,msize);
    TestMatrix *= 0.0;
    for (size_t i = 0; i < msize; ++i)
      {
        TestMatrix(i,i) = 1.0;
      }
    jiba::rvec s;
    jiba::rmat u,vt;
    jiba::SVD(TestMatrix,s,u,vt);
    for (size_t i = 0; i < msize; ++i)
      {
        BOOST_CHECK_CLOSE(s(i),1.0,std::numeric_limits<float>::epsilon());
        BOOST_CHECK_CLOSE(u(i,i),1.0,std::numeric_limits<float>::epsilon());
        BOOST_CHECK_CLOSE(vt(i,i),1.0,std::numeric_limits<float>::epsilon());
      }

  }

BOOST_AUTO_TEST_CASE(invert_matrix_test)
  {
    const size_t msize = 4;
    jiba::rmat TestMatrix(msize, msize);
    jiba::rmat InverseMatrix(msize, msize);
    srand(time(NULL));
    for (size_t i = 0; i < msize; ++i)
      for (size_t j = 0; j < msize; ++j)
        TestMatrix(i, j) = double(rand());

    jiba::InvertMatrix(TestMatrix, InverseMatrix);
    jiba::rmat MultMatrix(boost::numeric::ublas::prec_prod(TestMatrix,
        InverseMatrix));
    for (size_t i = 0; i < msize; ++i)
      {
        BOOST_CHECK_CLOSE(MultMatrix(i, i), 1.0,
            std::numeric_limits<float>::epsilon());
      }
  }

BOOST_AUTO_TEST_CASE(generelized_inverse_test)
  {
    const size_t msize = 4;
    jiba::rmat TestMatrix(msize, msize);
    jiba::rmat InverseMatrix(msize, msize);
    srand(time(NULL));
    for (size_t i = 0; i < msize; ++i)
      for (size_t j = 0; j < msize; ++j)
        TestMatrix(i, j) = double(rand());

    jiba::GeneralizedInverse(TestMatrix, InverseMatrix);
    jiba::rmat MultMatrix(boost::numeric::ublas::prec_prod(TestMatrix,
        InverseMatrix));
    for (size_t i = 0; i < msize; ++i)
      {
        BOOST_CHECK_CLOSE(MultMatrix(i, i), 1.0,
            std::numeric_limits<float>::epsilon());
      }
  }

BOOST_AUTO_TEST_SUITE_END()
