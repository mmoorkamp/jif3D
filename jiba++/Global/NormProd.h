//============================================================================
// Name        : NormProd.h
// Author      : May 30, 2010
// Version     : 
// Copyright   : 2010, mmoorkamp
//============================================================================


#ifndef NORMPROD_H_
#define NORMPROD_H_

#include <cassert>

/*! \file NormProd.h
 * Contains a template to calculate the inner product of two vectors given a norm by the diagonal of covariance matrix
 */

namespace jiba
  {
    /** \addtogroup util General utility routines */
    /* @{ */
    //! Calculate the inner product of two vectors given a norm by the diagonal of a covariance matrix
    /*! For inversion with respect to a model covariance we have to calculate the inner product
     * of two vectors with respect to the norm defined by that covariance. This function
     * implements this for diagonal covariances. These can be expressed in a vector as the third
     * argument. This can be easily implemented with ublas vector functions, however this function
     * is up to a factor two faster.
     * @param a The first vector for the inner product
     * @param b The second vector for the inner product
     * @param NormDiag The diagonal elements of the covariance matrix stored in a vector.
     * @return The inner product with respect to the norm
     */
    template<class UblasVec1, class UblasVec2, class UblasVec3>
    inline double NormProd(const UblasVec1 &a, const UblasVec2 &b,
        const UblasVec3 &NormDiag)
      {
        const size_t nelem = a.size();
        assert(nelem == b.size());
        assert(nelem == NormDiag.size());
        double result = 0.0;
        for (size_t i = 0; i < nelem; ++i)
          {
            result += a(i) * b(i) / NormDiag(i);
          }
        return result;
      }
  /* @} */
  }
#endif /* NORMPROD_H_ */
