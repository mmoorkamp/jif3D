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

    rmat RotateMTCovar(double angle, const rmat &InvCov, bool Inv)
      {
        using boost::numeric::ublas::prec_prod;
        //check that the matrix has the right size, 4x4
        const size_t nelem = 4;
        if (InvCov.size1() != nelem || InvCov.size2() != nelem)
          {
            throw FatalException(
                "Covariance matrix does not have correct size for rotation", __FILE__, __LINE__);
          }
        rmat Result(InvCov.size1(), InvCov.size2());
        //precompute the trigonometric values needed for the rotation
        const double cangle = cos(angle);
        const double sangle = sin(angle);
        //construct the matrix that transforms the covariances
        jif3D::rmat P(nelem, nelem);
        //the first row of the covariance matrix
        P(0, 0) = cangle * cangle;
        P(0, 1) = -cangle * sangle;
        P(0, 2) = P(0, 1);
        P(0, 3) = sangle * sangle;
        //all other rows are simply permutations (with sign changes)
        //of the first row
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

        //the covariance matrix and its inverse transform in a similar way
        // the inverse needs simply the transpose of the transformation matrix P
        // as P^T = P^-1
        if (Inv)
          {
            jif3D::rmat temp = prec_prod(InvCov, trans(P));
            Result = prec_prod(P, temp);
          }
        else
          {
            jif3D::rmat temp = prec_prod(InvCov, P);
            Result = prec_prod(trans(P), temp);
          }
        return Result;
      }

  }
