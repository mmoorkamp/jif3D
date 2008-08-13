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
    //! Calculate the singular value decomposition for a given matrix
    /*! Calculate the singular value decomposition using lapack
     * @param SensitivityMatrix The n*m input matrix, will be modified during the calculation
     * @param s The vector of singular values
     * @param u The n*n transformation matrix
     * @param vt The m*m transformation matrix
     */
    void SVD(rmat &SensitivityMatrix, rvec &s, rmat &u, rmat &vt)
      {
        const size_t size1 = SensitivityMatrix.size1();
        const size_t size2 = SensitivityMatrix.size2();
        s.resize(std::min(size1,size2));
        vt.resize(size2, size2);
        u.resize(size1, size1);
        boost::numeric::bindings::lapack::gesvd(SensitivityMatrix, s, u, vt);
      }

    bool InvertMatrix(const jiba::rmat& input, jiba::rmat& inverse);

    void GeneralizedInverse(const rmat &Input, rmat &Inverse, const double lowthresh = 0.0, const double upthresh = 1.0);

  }

#endif /* MATRIXTOOLS_H_ */
