//============================================================================
// Name        : test_rwcovar.cpp
// Author      : Max Moorkamp
// Version     :
// Copyright   : 2010, MM
//============================================================================

//test the file utility functions
#define BOOST_TEST_MODULE RWCovar test
#define BOOST_TEST_MAIN ...

#include <boost/test/included/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <cstdlib>
#include "ReadWriteCovariance.h"

BOOST_AUTO_TEST_SUITE( RWCovar_Test_Suite )

    BOOST_AUTO_TEST_CASE (test_readwrite)
      {
        srand(time(0));
        srand48(time(0));
        const size_t matsize = 10 + rand() % 50;
        jif3D::comp_mat CovMat(matsize, matsize);
        const size_t noffvalues = matsize + 2 * (rand() % 3 * matsize);
        for (size_t i = 0; i < matsize; ++i)
          {
            CovMat(i, i) = drand48();
          }
        for (size_t i = 0; i < noffvalues; ++i)
          {
            size_t index1 = rand() % matsize;
            size_t index2 = rand() % matsize;
            CovMat(index1, index2) = drand48();
            CovMat(index2, index1) = CovMat(index1, index2);
          }
        jif3D::WriteCovarianceToNetcdf("cov.nc", CovMat);

        jif3D::comp_mat Compare(matsize, matsize);
        jif3D::ReadCovarianceFromNetcdf("cov.nc", Compare);

        for (size_t i = 0; i < matsize; ++i)
          {
            for (size_t j = 0; j < matsize; ++j)
              {
                const double value1 = CovMat(i, j);
                const double value2 = Compare(i, j);
                BOOST_CHECK_CLOSE(value1, value2, 1e-5);
              }
          }

      }
    BOOST_AUTO_TEST_SUITE_END()
