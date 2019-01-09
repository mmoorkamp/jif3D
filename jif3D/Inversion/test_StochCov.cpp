//============================================================================
// Name        : test_optstep.cpp
// Author      : Apr 20, 2009
// Version     :
// Copyright   : 2009, mmoorkamp
//============================================================================

#define BOOST_TEST_MODULE StochCov test
#define BOOST_TEST_MAIN ...
#include "../Global/Jif3DTesting.h"
#include "../Global/VecMat.h"
#include <boost/test/included/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <stdlib.h>
#include <fstream>
#include "StochasticCovariance.h"
#include "DiagonalCovariance.h"
#include "../MT/X3DModel.h"

BOOST_AUTO_TEST_SUITE( StochCov_Test_Suite )

    BOOST_AUTO_TEST_CASE (StochCov_test)
      {
        jif3D::X3DModel Model;
        const size_t nx = 100;
        const size_t ny = 100;
        const size_t nz = 20;
        Model.SetMeshSize(nx, ny, nz);
        Model.SetHorizontalCellSize(1.0, 1.0, nx, ny);
        std::fill_n(Model.SetZCellSizes().begin(), nz, 1.0);

        double a = 5.0;
        double nu = 1.0;
        double sigma = 1.0;

        jif3D::StochasticCovariance Cov(nx,ny,nz, a, nu, sigma);
        jif3D::rvec m(nx * ny * nz, 0.0);
        std::generate(m.begin(), m.end(), drand48);
        jif3D::rvec mCm = Cov.ApplyCovar(m);

        jif3D::rvec mCmi = Cov.ApplyInvCovar(mCm);

        for (size_t i = 0; i < nx * ny * nz; ++i)
          if (std::abs(m(i)) > 1e-6)
            {
              //BOOST_CHECK_CLOSE(m(i), mCmi(i), 0.1);
            }

        const size_t ntries = 0;
        for (size_t i = 0; i < ntries; ++i)
          {
            std::generate(m.begin(), m.end(), drand48);
            jif3D::rvec mCm = Cov.ApplyCovar(m);
            jif3D::rvec mCmi = Cov.ApplyInvCovar(mCm);
          }
      }

    BOOST_AUTO_TEST_CASE (DiagCov_test)
      {
        const size_t nmod = 751;
        jif3D::rvec Covvec(nmod), m(nmod);
        for (double &c : Covvec)
          {
            c = std::abs(drand48()) + 0.01;
          }
        jif3D::DiagonalCovariance Cov(Covvec);

        std::generate(m.begin(),m.end(),drand48);
        jif3D::rvec cm = Cov.ApplyCovar(m);
        jif3D::rvec cmi = Cov.ApplyInvCovar(cm);
        for (size_t i = 0; i < nmod; ++i)
          {
            BOOST_CHECK_CLOSE(m(i), cmi(i), 0.1);
          }
      }

    BOOST_AUTO_TEST_SUITE_END()
