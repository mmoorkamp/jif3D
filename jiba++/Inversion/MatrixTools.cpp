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
    namespace ublas = boost::numeric::ublas;
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
        G.resize(Input.size1(),Input.size2());
        G.assign(Input);
        //do the SVD
        SVD(G, s, u, vt);
        //translate the relative threshold values into absolute numbers
        const double maxeval = *std::max_element(s.begin(), s.end());
        const double absupper = maxeval * upthresh;
        const double abslower = maxeval * lowthresh;
        const size_t neval = s.size();
        //this will eventually hold lambda^-1 * vt
        rmat lambdavt(vt.size1(),neval);
        // go through all eigenvalues
        for (size_t i = 0; i < neval; ++i)
          {
            //check if we need to apply the threshold
            if (s(i)> absupper || s(i) < abslower)
              {
                //if so make the eigenvectors zero
                column(u, i) *= 0.0;
                std::fill_n(column(lambdavt, i).begin(),vt.size1(), 0.0);
              }
            else
              {
                //otherwise devide vT by the eigenvalue
                // this is the multiplication lambda^-1 vT
                column(lambdavt, i) = 1. / s(i) * row(vt,i);
              }
          }
        //alocate space for the inverse
        Inverse.resize(Input.size2(), Input.size1());
        //do the multiplication (v lambda^-1) * uT
        atlas::gemm(CblasNoTrans, CblasTrans, 1.0, lambdavt, u, 0.0, Inverse);
      }

    /*!
     * Returns the determinant of a real square Matrix. This is
     * a general method that uses LU-factorization and therefore
     * works for any square matrix.
     * @param Matrix The real square matrix
     * @return det(Matrix)
     */
    double Determinant(const jiba::rmat &Matrix)
      {
        assert(Matrix.size1() == Matrix.size2());
        jiba::rmat Factor(Matrix);
        typedef ublas::permutation_matrix<std::size_t> pmatrix;
        pmatrix pm(Factor.size1());
        lu_factorize(Factor,pm);

        double det = 1.0;
        for (size_t i = 0; i < Factor.size1(); ++i)
          {
            if (pm(i) != i)
              det *= -Factor(i,i);
            else
              det *= Factor(i,i);
          }
        return det;
      }
  }
