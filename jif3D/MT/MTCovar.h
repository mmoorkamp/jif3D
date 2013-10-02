//============================================================================
// Name        : MTCovar.h
// Author      : 11 Jul 2013
// Version     : 
// Copyright   : 2013, mm489
//============================================================================

#ifndef MTCOVAR_H_
#define MTCOVAR_H_

#include <boost/numeric/ublas/symmetric.hpp>
#include "../Global/VecMat.h"

namespace jif3D
  {
    /** \addtogroup mtmodelling Forward modelling of magnetotelluric data */
    /* @{ */
    typedef boost::numeric::ublas::symmetric_matrix<double> covmat;

    rmat RotateMTCovar(double angle,const rmat &Cov, bool Inv = false);
  /* @} */
  }

#endif /* MTCOVAR_H_ */
