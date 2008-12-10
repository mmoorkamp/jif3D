//============================================================================
// Name        : BasicGravElements.h
// Author      : Dec 10, 2008
// Version     :
// Copyright   : 2008, mmoorkamp
//============================================================================


#ifndef BASICGRAVELEMENTS_H_
#define BASICGRAVELEMENTS_H_

#include "ThreeDGravityModel.h"

namespace jiba
  {
    //! The constant of gravity
    static const double Grav_const = 6.67428e-8; // in units cm^3/g s

    //! Calculate a single geometric term for the graviational acceleration due to a rectangular prism
    double CalcGravTerm(const double x, const double y, const double z);
    //! Calculate the geometric term  due to a rectangular prism
    double CalcGravBoxTerm(const double meas_x, const double meas_y,
        const double meas_z, const double ul_corner_x,
        const double ul_corner_y, const double ul_corner_z,
        const double x_size, const double y_size, const double z_size);

    //! Calculate the geometric part of the gravimetry matrix for a single rectangular prism
    GravimetryMatrix CalcTensorBoxTerm(const double meas_x,
        const double meas_y, const double meas_z, const double ul_corner_x,
        const double ul_corner_y, const double ul_corner_z,
        const double x_size, const double y_size, const double z_size);

    double CalcUxxTerm(const double meas_x, const double meas_y,
        const double meas_z, const double ul_corner_x,
        const double ul_corner_y, const double ul_corner_z,
        const double x_size, const double y_size, const double z_size);

    double CalcUxyTerm(const double meas_x, const double meas_y,
        const double meas_z, const double ul_corner_x,
        const double ul_corner_y, const double ul_corner_z,
        const double x_size, const double y_size, const double z_size);

    //! the gravitational acceleration of an semi-infinite slab
    double CalcGravSemiInfSheet(const double hor_dist, const double ver_dist,
        const double thick, const double density);
    //! Calculate one of the terms of diagonal elements of the gravimetric matrxi
    double CalcFTGDiagonalTerm(const double a, const double b, const double c);
    //! Calculate one of the terms of off-diagonal elements of the gravimetric matrxi
    double CalcFTGOffDiagonalTerm(const double value, const double x,
        const double y, const double z);
    //! the "geometric term" for gravitational acceleration of an infinite slab
    inline double CalcInfSheetTerm(const double thick)
      {
        return 2.0 * M_PI * Grav_const * thick;
      }
  //end of namespace jiba
  }
#endif /* BASICGRAVELEMENTS_H_ */
