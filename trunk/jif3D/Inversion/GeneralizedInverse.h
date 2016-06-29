//============================================================================
// Name        : GeneralizedInverse.h
// Author      : 17 May 2012
// Version     :
// Copyright   : 2012, mm489
//============================================================================

#ifndef GENERALIZEDINVERSE_H_
#define GENERALIZEDINVERSE_H_

#include "../Global/VecMat.h"
#include "../Global/Jif3DGlobal.h"

/*! \file GeneralizedInverse.h
 * This file includes  operations on matrices that are associated
 * with singular value decomposition and generalized inverse.
 * The operations are implemented in terms of calls to blas, atlas through
 * the ublas numeric bindings and therefore introduce a number of dependencies.
 */

namespace jif3D
  {

    /** \addtogroup inversion General routines for inversion */
    /* @{ */

    //! Calculate the singular value decomposition for a given matrix
    /*! Calculate the singular value decomposition \f$ M = uSv^T \f$ using lapack. As \f$ S\f$ only has non-zero
     * elements on its  diagonal, it is stored as a vector.
     * @param SensitivityMatrix The n*m input matrix, will be modified during the calculation
     * @param s The vector of singular values
     * @param u The n*n transformation matrix
     * @param vt The m*m transformation matrix
     */
    J3DEXPORT void SVD(rmat &SensitivityMatrix, rvec &s, rmat &u, rmat &vt);

    //! Implements the generalized inverse through truncated SVD
    /*! This class provides the functionality to calculate a
     * generalized inverse based on the truncated SVD and related
     * quantities such as model resolution and covariance.
     */
    class J3DEXPORT GeneralizedInverse
      {
    private:
      //! The singular values of the matrix G
      rvec s;
      //! The data space eigenvectors of the matrix G
      rmat u;
      //! The model space eigenvectors of the matrix G
      rmat vt;
      //! The input matrix to the SVD gets overwritten by the algorithm, as we do not want this we store a copy of the sensitivity matrix
      rmat G;
    public:
      //! return the vector of generalized eigenvalues, no threshold is applied
      const rvec &GetEigenvalues() const
        {
          return s;
        }
      //! return the data eigenvectors after application of the thresholds
      const rmat &GetDataVectors() const
        {
          return u;
        }
      //! return the model eigenvectors after application of the thresholds
      const rmat &GetModelVectors() const
        {
          return vt;
        }
      //! Return the resolution matrix after application of the thresholds
      rmat GetModelResolution() const
        {
          return ublas::prec_prod(ublas::trans(vt), vt);
        }
      //rmat GetModelCovariance() const {return ublas::prec_prod(SensitivityInv,ublas::trans(SensitivityInv));}
      //! The core calculation routine to calculate the generalized inverse
      void operator()(const rmat &Input, rmat &Inverse, const double lowthresh = 0.0,
          const double upthresh = 1.0);
      //! The constructor just initialises everything to the default
      GeneralizedInverse() :
          s(), u(), vt(), G()
        {

        }
      };
  /* @} */
  }
#endif /* GENERALIZEDINVERSE_H_ */
