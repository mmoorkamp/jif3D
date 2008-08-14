//============================================================================
// Name        : MatrixTools.cpp
// Author      : Max Moorkamp
// Version     :
// Copyright   : 2008, MM
//============================================================================


#include "MatrixTools.h"
#include <boost/numeric/bindings/atlas/cblas3.hpp>

namespace jiba
  {
    namespace atlas = boost::numeric::bindings::atlas;

    void GeneralizedInverse(const rmat &Input, rmat &Inverse,
        const double lowthresh, const double upthresh)
      {
        rvec s;
        rmat u, vt;
        rmat G(Input);
        SVD(G, s, u, vt);
        const double maxeval = *std::max_element(s.begin(), s.end());
        const double absupper = maxeval * upthresh;
        const double abslower = maxeval * lowthresh;
        const size_t neval = s.size();
        for (size_t i = 0; i < neval; ++i)
          {
            if (s(i) > absupper || s(i) < abslower)
              {
                column(u, i) *= 0.0;
                row(vt, i) *= 0.0;
              }
            else
              {
                row(vt, i) *= 1. / s(i);
              }
          }

        boost::numeric::ublas::matrix_range<jiba::rmat> nonzero(vt,
            boost::numeric::ublas::range(0, neval),
            boost::numeric::ublas::range(0, vt.size1()));

        Inverse.resize(Input.size2(), Input.size1());
        atlas::gemm(CblasTrans, CblasTrans, 1.0, nonzero, u, 0.0, Inverse);
      }

    bool InvertMatrix(const jiba::rmat& input, jiba::rmat& inverse)
      {
        using namespace boost::numeric::ublas;
        typedef permutation_matrix<std::size_t> pmatrix;
        // create a working copy of the input
        jiba::rmat A(input);
        // create a permutation matrix for the LU-factorization
        pmatrix pm(A.size1());

        // perform LU-factorization
        int res = lu_factorize(A, pm);
        if (res != 0)
          return false;

        // create identity matrix of "inverse"
        inverse.assign(boost::numeric::ublas::identity_matrix<double>(A.size1()));

        // backsubstitute to get the inverse
        lu_substitute(A, pm, inverse);

        return true;
      }
  }
