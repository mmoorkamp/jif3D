//============================================================================
// Name        : BasicGravElements.h
// Author      : Dec 10, 2008
// Version     :
// Copyright   : 2008, mmoorkamp
//============================================================================

#ifndef BASICGRAVELEMENTS_H_
#define BASICGRAVELEMENTS_H_


#include "../Global/VecMat.h"
#include "../Global/NumUtil.h"
#include "../Global/Jif3DGlobal.h"
#include <boost/math/constants/constants.hpp>

namespace jif3D
  {
    /** \addtogroup gravity Gravity forward modeling, display and inversion */
    /* @{ */
    //! The constant of gravity  in units m^3/kg s
    static const double Grav_const = 6.67428e-11; // in units m^3/kg s

    //! Calculate the geometric term for scalar gravity for a rectangular prism
    J3DEXPORT double CalcGravBoxTerm(const double meas_x, const double meas_y,
        const double meas_z, const double ul_corner_x,
        const double ul_corner_y, const double ul_corner_z,
        const double x_size, const double y_size, const double z_size);

    //! Calculate the geometric part of the gravimetry matrix for a single rectangular prism
    J3DEXPORT rmat CalcTensorBoxTerm(const double meas_x,
        const double meas_y, const double meas_z, const double ul_corner_x,
        const double ul_corner_y, const double ul_corner_z,
        const double x_size, const double y_size, const double z_size);
    //! Calculate the Uxx geometric term for a single rectangular prism, is also used for the other diagonal elements with permuted arguments
    J3DEXPORT double CalcUxxTerm(const double meas_x, const double meas_y,
        const double meas_z, const double ul_corner_x,
        const double ul_corner_y, const double ul_corner_z,
        const double x_size, const double y_size, const double z_size);
    //! Calculate the Uxy geometric term for a single rectangular prism, is also used for the other off-diagonal elements with permuted arguments
    J3DEXPORT double CalcUxyTerm(const double meas_x, const double meas_y,
        const double meas_z, const double ul_corner_x,
        const double ul_corner_y, const double ul_corner_z,
        const double x_size, const double y_size, const double z_size);

    //! Calculate the scalar geometric term for gravitational acceleration of an semi-infinite slab
    J3DEXPORT double CalcGravSemiInfSheet(const double hor_dist, const double ver_dist,
        const double thick, const double density);

    //! Calculate the scalar geometric term for gravitational acceleration of an infinite slab
    inline double CalcInfSheetTerm(const double measz, const double top,
        const double bottom)
      {
        const double thick = bottom - top;
        //we have to distinguish two cases
        //the measurement is outside the layer
        if ((measz <= top) || (measz >= bottom))
          {
            return 2.0 * boost::math::constants::pi<double>() * Grav_const * thick * jif3D::sign(top - measz);
          }
        //if we get here it means the measurement is inside the layer
        return 2.0 * boost::math::constants::pi<double>() * Grav_const * (thick - 2.0 * (measz - top));
      }
    //! An infinite sheet has only an effect on \f$ U_{zz}\f$ if the measurement is take within the boundaries of that layer
    inline double CalcUzzInfSheetTerm(const double measz, const double top,
        const double bottom)
      {
        //outside the layer the effect is zero
        if ((measz <= top) || (measz >= bottom))
          {
            return 0.0;
          }
        //if we get here we are inside the layer
        return 4.0 * boost::math::constants::pi<double>() * Grav_const;
      }
  /* @} */
  //end of namespace jif3D
  }
#endif /* BASICGRAVELEMENTS_H_ */
