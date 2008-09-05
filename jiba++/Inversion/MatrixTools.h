//============================================================================
// Name        : MatrixTools.h
// Author      : Max Moorkamp
// Version     :
// Copyright   : 2008, MM
//============================================================================

#ifndef MATRIXTOOLS_H_
#define MATRIXTOOLS_H_

#include "../Global/VecMat.h"


namespace jiba
  {
    namespace ublas = boost::numeric::ublas;
    /** \addtogroup inversion General routines for inversion */
    /* @{ */

    //! Calculate the singular value decomposition for a given matrix
    /*! Calculate the singular value decomposition using lapack
     * @param SensitivityMatrix The n*m input matrix, will be modified during the calculation
     * @param s The vector of singular values
     * @param u The n*n transformation matrix
     * @param vt The m*m transformation matrix
     */
    inline void SVD(rmat &SensitivityMatrix, rvec &s, rmat &u, rmat &vt)
      {
        const size_t size1 = SensitivityMatrix.size1();
        const size_t size2 = SensitivityMatrix.size2();
        s.resize(std::min(size1, size2));
        vt.resize(size2, size2);
        u.resize(size1, size1);
        boost::numeric::bindings::lapack::gesvd(SensitivityMatrix, s, u, vt);
      }
    //! Invert a square Matrix using LU-Factorization
    bool InvertMatrix(const jiba::rmat& input, jiba::rmat& inverse);
    //! Calculate the determinant of a real square matrix
    double Determinant(const jiba::rmat &Matrix);
    //! Implements the generalized inverse through truncated SVD
    class GeneralizedInverse
      {
    private:
      rvec s;
      rmat u;
      rmat vt;
      rmat G;
    public:
      //! return the vector of generalized eigenvalues, no threshold is applied
      const rvec &GetEigenvalues() const {return s;}
      //! return the data eigenvectors after application of the thresholds
      const rmat &GetDataVectors() const {return u;}
      //! return the model eigenvectors after application of the thresholds
      const rmat &GetModelVectors() const {return vt;}
      //! Return the resolution matrix after application of the thresholds
      rmat GetModelResolution() const {return ublas::prec_prod(ublas::trans(vt),vt);}
      //rmat GetModelCovariance() const {return ublas::prec_prod(SensitivityInv,ublas::trans(SensitivityInv));}
      //! The core calculation routine to calculate the generalized inverse
      void operator()(const rmat &Input, rmat &Inverse, const double lowthresh = 0.0,
          const double upthresh = 1.0);
      };
    /* @} */
  }

#endif /* MATRIXTOOLS_H_ */
