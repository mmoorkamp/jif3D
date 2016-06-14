//============================================================================
// Name        : MatrixTools.cpp
// Author      : Max Moorkamp
// Version     :
// Copyright   : 2008, MM
//============================================================================


#include "MatrixTools.h"

namespace jif3D
  {

    namespace ublas = boost::numeric::ublas;

    /*! Returns the determinant of a real square Matrix. This is
     * a general method that uses LU-factorization and therefore
     * works for any square matrix.
     * @param Matrix The real square matrix
     * @return det(Matrix)
     */
    double Determinant(const jif3D::rmat &Matrix)
      {
        assert(Matrix.size1() == Matrix.size2());
        jif3D::rmat Factor(Matrix);
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
