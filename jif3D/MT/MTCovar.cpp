//============================================================================
// Name        : MTCovar.cpp
// Author      : 11 Jul 2013
// Version     : 
// Copyright   : 2013, mm489
//============================================================================

#include "MTCovar.h"
#include "../Global/FatalException.h"
namespace jif3D
  {

    covmat RotateMTCovar(double angle, const covmat &Cov)
      {
        const size_t nelem = 4;
        if (Cov.size1() != nelem || Cov.size2() != nelem)
          {
            throw FatalException(
                "Covariance matrix does not have correct size for rotation");
          }
        covmat Result(Cov.size1(), Cov.size2());
        const double cangle = cos(angle);
        const double sangle = sin(angle);

        jif3D::rmat P(Cov.size1(), Cov.size2());
        P(0, 0) = cangle * cangle;
        P(0, 1) = cangle * sangle;
        P(0, 2) = P(0, 1);
        P(0, 3) = sangle * sangle;
        P(1, 0) = -P(0, 1);
        P(1, 1) = P(0, 0);
        P(1, 2) = -P(0, 3);
        P(1, 3) = P(0, 1);
        P(2, 0) = -P(0, 1);
        P(2, 1) = -P(0, 3);
        P(2, 2) = P(0, 0);
        P(2, 3) = P(0, 1);
        P(3, 0) = P(0, 3);
        P(3, 1) = -P(0, 1);
        P(3, 2) = -P(0, 1);
        P(3, 3) = P(0, 0);
        jif3D::rmat temp = boost::numeric::ublas::prec_prod(Cov, P);
        Result = boost::numeric::ublas::prec_prod(trans(P), temp);
        return Result;
      }

  }
