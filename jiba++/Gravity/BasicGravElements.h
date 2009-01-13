//============================================================================
// Name        : BasicGravElements.h
// Author      : Dec 10, 2008
// Version     :
// Copyright   : 2008, mmoorkamp
//============================================================================


#ifndef BASICGRAVELEMENTS_H_
#define BASICGRAVELEMENTS_H_

#include "ThreeDGravityModel.h"
#include "../Global/NumUtil.h"
namespace jiba
  {
    /** \addtogroup gravity Gravity forward modelling, display and inversion */
    /* @{ */
    //! The constant of gravity
    static const double Grav_const = 6.67428e-8; // in units cm^3/g s

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
    //! Calculate the Uxx geometric term for a single rectangular prism, is also used for the other diagonal elements with permuted arguments
    double CalcUxxTerm(const double meas_x, const double meas_y,
        const double meas_z, const double ul_corner_x,
        const double ul_corner_y, const double ul_corner_z,
        const double x_size, const double y_size, const double z_size);
    //! Calculate the Uxy geometric term for a single rectangular prism, is also used for the other off-diagonal elements with permuted arguments
    double CalcUxyTerm(const double meas_x, const double meas_y,
        const double meas_z, const double ul_corner_x,
        const double ul_corner_y, const double ul_corner_z,
        const double x_size, const double y_size, const double z_size);

    //! Calculate the "geometric term" for gravitational acceleration of an semi-infinite slab
    double CalcGravSemiInfSheet(const double hor_dist, const double ver_dist,
        const double thick, const double density);

    //! Calculate the "geometric term" for gravitational acceleration of an infinite slab
    inline double CalcInfSheetTerm(const double measz, const double top,
        const double bottom)
      {
        const double thick = bottom - top;
        if ((measz <= top) || (measz >= bottom))
          {
            return 2.0 * M_PI * Grav_const * thick * jiba::sign(top - measz);
          }
        return 2.0 * M_PI * Grav_const * (thick - 2.0 * (measz - top));
      }

    inline double CalcUzzInfSheetTerm(const double measz, const double top,
        const double bottom)
      {
        if ((measz <= top) || (measz >= bottom))
          {
            return 0.0;
          }
        return 4.0 * M_PI * Grav_const;
      }
  /* @} */
  //end of namespace jiba
  }
#endif /* BASICGRAVELEMENTS_H_ */
