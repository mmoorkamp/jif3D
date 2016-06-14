//============================================================================
// Name        : MTCovar.h
// Author      : 11 Jul 2013
// Version     : 
// Copyright   : 2013, mm489
//============================================================================

#ifndef MTCOVAR_H_
#define MTCOVAR_H_


#include "../Global/VecMat.h"
#include <boost/numeric/ublas/symmetric.hpp>

namespace jif3D
  {
    /** \addtogroup mtmodelling Forward modelling of magnetotelluric data */
    /* @{ */
    typedef boost::numeric::ublas::symmetric_matrix<double> covmat;

    //! Rotate a covariance matrix containing error information for MT data
    /*! When rotating MT impedance data from one coordinate system to another, we
     * have to consider not only the standard error, i.e. the variances, but the
     * full covariance matrix. This function performs the associated rotation of
     * a covariance matrix from one coordinate system to another. As the rotation
     * acts on each impedance tensor separately we do not rotate the martrix for
     * all impedances at all frequencies and stations, but just the 4x4 matrix
     * for a single tensor.
     * @param angle The rotation angle in radian
     * @param Cov The 4x4 matrix containing the covariance information
     * @param Inv Do we rotate the covariance matrix or its inverse, the weighting matrix
     * @return 4x4 covariance matrix in the rotated coordianate system
     */
    rmat RotateMTCovar(double angle,const rmat &Cov, bool Inv = false);
  /* @} */
  }

#endif /* MTCOVAR_H_ */
