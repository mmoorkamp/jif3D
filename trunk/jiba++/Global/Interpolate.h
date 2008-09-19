//============================================================================
// Name        : Interpolate.h
// Author      : Max Moorkamp
// Version     :
// Copyright   : 2008, MM
//============================================================================

#ifndef INTERPOLATE_H_
#define INTERPOLATE_H_

#include "VecMat.h"

namespace jiba
  {

    /** \addtogroup util General utility routines */
    /* @{ */

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
