//============================================================================
// Name        : MatrixTools.h
// Author      : Max Moorkamp
// Version     :
// Copyright   : 2008, MM
//============================================================================

#ifndef MATRIXTOOLS_H_
#define MATRIXTOOLS_H_

#include "../Global/VecMat.h"

/*! \file MatrixTools.h
 * This file includes a few mathematical operations on matrices
 * like  matrix inversion and calculation of the determinant.
 */

namespace jif3D
  {
    namespace ublas = boost::numeric::ublas;
    /** \addtogroup inversion General routines for inversion */
    /* @{ */

    //! Invert a square Matrix using LU-Factorization
    /*! This function inverts a square real matrix using LU-factorization
     * @param input The original matrix
     * @param inverse The inverse
     * @return True if success, false otherwise
     */
    template<class InMatrix, class OutMatrix>
    bool InvertMatrix(InMatrix& input, OutMatrix& inverse)
      {
        using namespace boost::numeric::ublas;
        typedef permutation_matrix<std::size_t> pmatrix;

        // create a permutation matrix for the LU-factorization
        pmatrix pm(input.size1());

        // perform LU-factorization
        int res = lu_factorize(input, pm);
        if (res != 0)
          return false;

        // create identity matrix of "inverse"
        inverse.assign(identity_matrix<double>(input.size1()));

        // backsubstitute to get the inverse
        lu_substitute(input, pm, inverse);

        return true;
      }
    //! Calculate the determinant of a real square matrix
    double Determinant(const jif3D::rmat &Matrix);

  /* @} */
  }

#endif /* MATRIXTOOLS_H_ */
