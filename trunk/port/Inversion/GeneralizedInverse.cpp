//============================================================================
// Name        : GeneralizedInverse.cpp
// Author      : 17 May 2012
// Version     : 
// Copyright   : 2012, mm489
//============================================================================

#include "GeneralizedInverse.h"
#include <boost/numeric/bindings/lapack/lapack.hpp>
#include <boost/numeric/bindings/lapack/geev.hpp>
#include <boost/numeric/bindings/traits/ublas_matrix.hpp>
#include <boost/numeric/bindings/traits/ublas_vector.hpp>
#include <boost/numeric/bindings/traits/vector_traits.hpp>
#include <boost/numeric/bindings/atlas/cblas3.hpp>

namespace jif3D
  {
    namespace atlas = boost::numeric::bindings::atlas;
    namespace ublas = boost::numeric::ublas;

    void SVD(rmat &SensitivityMatrix, rvec &s, rmat &u, rmat &vt)
      {
        const size_t size1 = SensitivityMatrix.size1();
        const size_t size2 = SensitivityMatrix.size2();
        s.resize(std::min(size1, size2));
        vt.resize(size2, size2);
        u.resize(size1, size1);
        boost::numeric::bindings::lapack::gesvd(SensitivityMatrix, s, u, vt);
      }

    /*! The routine that calculates the generalized inverse for a given matrix
     * @param Input The \f$n\times m\f$ matrix we want to invert
     * @param Inverse The \f$ m \times n \f$ result of the inversion
     * @param lowthresh The lower truncation threshold relative to the maximum eigenvalue
     * @param upthresh The upper truncation threshold relative to the maximum eigenvalue, should usually be 1
     */
    void GeneralizedInverse::operator()(const rmat &Input, rmat &Inverse,
        const double lowthresh, const double upthresh)
      {
        //allocate memory
        G.resize(Input.size1(), Input.size2());
        G.assign(Input);
        //do the SVD
        SVD(G, s, u, vt);
        //translate the relative threshold values into absolute numbers
        const double maxeval = *std::max_element(s.begin(), s.end());
        const double absupper = maxeval * upthresh;
        const double abslower = maxeval * lowthresh;
        const size_t neval = s.size();
        //this will eventually hold lambda^-1 * vt
        rmat lambdavt(vt.size1(), neval);
        // go through all eigenvalues
        for (size_t i = 0; i < neval; ++i)
          {
            //check if we need to apply the threshold
            if (s(i) > absupper || s(i) < abslower)
              {
                //if so make the eigenvectors zero
                column(u, i) *= 0.0;
                std::fill_n(column(lambdavt, i).begin(), vt.size1(), 0.0);
              }
            else
              {
                //otherwise devide vT by the eigenvalue
                // this is the multiplication lambda^-1 vT
                column(lambdavt, i) = 1. / s(i) * row(vt, i);
              }
          }
        //alocate space for the inverse
        Inverse.resize(Input.size2(), Input.size1());
        //do the multiplication (v lambda^-1) * uT
        atlas::gemm(CblasNoTrans, CblasTrans, 1.0, lambdavt, u, 0.0, Inverse);
      }
  }

