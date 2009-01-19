//============================================================================
// Name        : Interpolate.h
// Author      : Max Moorkamp
// Version     :
// Copyright   : 2008, MM
//============================================================================

#ifndef INTERPOLATE_H_
#define INTERPOLATE_H_

#include "VecMat.h"

/*! \file Interpolate.h
 * This file contains routines for interpolation and extrapolation.
 */
namespace jiba
  {

    /** \addtogroup util General utility routines */
    /* @{ */
    /*! Perform quadratic interpolation for the minimum given three
     * pairs of function values \f$ x_i, y_i=f(x_i)\f$
     * @param x1 The first x-value
     * @param y1 The associated function value f(x1)
     * @param x2 The second x-value
     * @param y2 The associated function value f(x2)
     * @param x3 The third x-value
     * @param y3 The associated function value f(x3)
     * @return The x-value of the interpolated minimum
     */
    double QuadraticInterpolation(const double x1, const double y1,
        const double x2, const double y2, const double x3, const double y3)
      {
        double value = (x2 * x2 - x3 * x3) * y1 + (x3 * x3 - x1 * x1) * y2
            + (x1 * x1 - x2 * x2) * y3;
        value /= (x2 - x3) * y1 + (x3 - x1) * y2 + (x1 - x2) * y3;
        value /= 2.0;
        return value;

      }
  /* @} */
  }
#endif /* INTERPOLATE_H_ */
