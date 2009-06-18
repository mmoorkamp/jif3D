//============================================================================
// Name        : BasicGravElements.cpp
// Author      : Dec 10, 2008
// Version     :
// Copyright   : 2008, mmoorkamp
//============================================================================

#include <cmath>
#include <algorithm>
#include <numeric>
#include <limits>
#include "BasicGravElements.h"

namespace jiba
  {
    //! Calculate a single geometric term for the gravitational acceleration due to a rectangular prism
    /*! Calculate one term for the gravitational potential of a box, we use the nomenclature of eq. 4-6 in Li and Chouteau.
     * The parameters x,y and z are the distances to the corners of the box, we will call this functions with different
     * permutations of real world x,y and z coordinates, so only in some cases x corresponds to the x-axis.
     * The various terms all have the form \f$ x * log(y +r) + y * log(x+r) + z*atan2(z*r, x*y) \f$
     * where x,y and z are distances in x, y and z-direction for the current corner, respectively and
     * \f$ r = \sqrt{x^2 +y^2 +z^2} \f$.
     */
    inline double CalcGravTerm(const double x, const double y, const double z)
      {
        //if the distance r between the measurement point and one of the border points is very small
        //the log(...) expressions become undefined, so we shift the measurement point by a tiny amount
        const double r = std::max(sqrt(x * x + y * y + z * z),
            std::numeric_limits<double>::epsilon());
        double rvalue = x * log(y + r) + y * log(x + r);
        // atan2 takes care of small denominators
        rvalue += z * atan2(z * r, x * y);

        return rvalue;
      }

    //! Calculate one of the terms in the calculation of diagonal elements of the gravimetric matrxi
    /*! The terms \f$ U_{xx}, U_{yy} \f$ and \f$ U_{zz} \f$  of the gravimetric matrix
     * all are sums of terms of the form \f$ atan \frac{a *b}{c *r } \f$
     */
    inline double CalcFTGDiagonalTerm(const double a, const double b,
        const double c)
      {
        return atan2(a * b, c * sqrt(a * a + b * b + c * c));
      }

    //! Calculate one of the terms of off-diagonal elements of the gravimetric matrixs
    /*! The terms \f$ U_{xy}, U_{xz} \f$ and \f$ U_{yz} \f$  of the gravimetric matrix
     * all are sums of terms of the form \f$ \log (x +r) \f$
     */
    inline double CalcFTGOffDiagonalTerm(const double value, const double x,
        const double y, const double z)
      {
        //if the logarithmic term becomes really small we cheat
        return log(std::max(value + sqrt(x * x + y * y + z * z),
            std::numeric_limits<double>::epsilon()));
      }

    /*! Calculate the geometric term for the gravitational potential of a rectangular prism at a point meas_{x,y,z}
     *  The calculation is based on equation 4 in Li and Chouteau, Surveys in Geophysics, 19, 339-368, 1998
     * Given the coordinates of the measurements (meas_x, meas_y and meas_z), the
     * coordinates of the upper front left corner (ul_corner_x etc.), the size
     * in the three directions (x_size, y_size and z_size) and the density, we calculate
     * the geometric part caused by a prism with these parameters. All dimensions are in meters.
     * If we multiply the result of this function with density
     * in \f$ g/cm^3 \f$ we get the acceleration in \f$m/s^2\f$.
     *
     * This equation works as long as the measurement point is not on one of the corners of the box.
     * @param meas_x x-coordinate of the measurement in m
     * @param meas_y y-coordinate of the measurement in m
     * @param meas_z z-coordinate of the measurement in m
     * @param ul_corner_x x-coordinate of the upper left front corner of the density element in m
     * @param ul_corner_y y-coordinate of the upper left front corner of the density element in m
     * @param ul_corner_z z-coordinate of the upper left front corner of the density element in m
     * @param x_size size of the prism in x-direction in m
     * @param y_size size of the prism in y-direction in m
     * @param z_size size of the prism in z-direction in m
     * @return Gravitational acceleration in m/s^2
     */
    double CalcGravBoxTerm(const double meas_x, const double meas_y,
        const double meas_z, const double ul_corner_x,
        const double ul_corner_y, const double ul_corner_z,
        const double x_size, const double y_size, const double z_size)
      {
        //we need 8 terms, one for each corner of the prism, they have varying signs
        double returnvalue = -CalcGravTerm(meas_x - ul_corner_x, meas_y
            - ul_corner_y, meas_z - ul_corner_z);
        // one corner shifted -> positive sign
        returnvalue += CalcGravTerm(meas_x - (ul_corner_x + x_size), meas_y
            - ul_corner_y, meas_z - ul_corner_z);
        returnvalue += CalcGravTerm(meas_x - ul_corner_x, meas_y - (ul_corner_y
            + y_size), meas_z - ul_corner_z);
        returnvalue += CalcGravTerm(meas_x - ul_corner_x, meas_y - ul_corner_y,
            meas_z - (ul_corner_z + z_size));
        // two corners shifted -> negative sign
        returnvalue -= CalcGravTerm(meas_x - (ul_corner_x + x_size), meas_y
            - (ul_corner_y + y_size), meas_z - ul_corner_z);
        returnvalue -= CalcGravTerm(meas_x - (ul_corner_x + x_size), meas_y
            - ul_corner_y, meas_z - (ul_corner_z + z_size));
        returnvalue -= CalcGravTerm(meas_x - ul_corner_x, meas_y - (ul_corner_y
            + y_size), meas_z - (ul_corner_z + z_size));
        //three corners shifted -> positive sign
        returnvalue += CalcGravTerm(meas_x - (ul_corner_x + x_size), meas_y
            - (ul_corner_y + y_size), meas_z - (ul_corner_z + z_size));
        // we multiply the geometric term by the gravitational constant and the density
        return -returnvalue * Grav_const;
      }

    /*! Calculate the geometric term of the xx-element of the second order derivative tensor, we use the same
     * function with permutated parameters to calculate the other two diagonal elements.
     * The description of the parameters is for the Uxx calculation.
     * @param meas_x x-coordinate of the measurement in m
     * @param meas_y y-coordinate of the measurement in m
     * @param meas_z z-coordinate of the measurement in m
     * @param ul_corner_x x-coordinate of the upper left front corner of the density element in m
     * @param ul_corner_y y-coordinate of the upper left front corner of the density element in m
     * @param ul_corner_z z-coordinate of the upper left front corner of the density element in m
     * @param x_size size of the prism in x-direction in m
     * @param y_size size of the prism in y-direction in m
     * @param z_size size of the prism in z-direction in m
     * @return The Uxx element of the FTG tensor
     */
    double CalcUxxTerm(const double meas_x, const double meas_y,
        const double meas_z, const double ul_corner_x,
        const double ul_corner_y, const double ul_corner_z,
        const double x_size, const double y_size, const double z_size)
      {
        //we need 8 terms, one for each corner of the prism, they have varying signs
        double returnvalue = -CalcFTGDiagonalTerm(meas_y - ul_corner_y, meas_z
            - ul_corner_z, meas_x - ul_corner_x);
        // one corner shifted -> positive sign
        returnvalue += CalcFTGDiagonalTerm(meas_y - ul_corner_y, meas_z
            - ul_corner_z, meas_x - (ul_corner_x + x_size));
        returnvalue += CalcFTGDiagonalTerm(meas_y - (ul_corner_y + y_size),
            meas_z - ul_corner_z, meas_x - ul_corner_x);
        returnvalue += CalcFTGDiagonalTerm(meas_y - ul_corner_y, meas_z
            - (ul_corner_z + z_size), meas_x - ul_corner_x);
        // two corners shifted -> negative sign
        returnvalue -= CalcFTGDiagonalTerm(meas_y - (ul_corner_y + y_size),
            meas_z - ul_corner_z, meas_x - (ul_corner_x + x_size));
        returnvalue -= CalcFTGDiagonalTerm(meas_y - ul_corner_y, meas_z
            - (ul_corner_z + z_size), meas_x - (ul_corner_x + x_size));
        returnvalue -= CalcFTGDiagonalTerm(meas_y - (ul_corner_y + y_size),
            meas_z - (ul_corner_z + z_size), meas_x - ul_corner_x);
        //three corners shifted -> positive sign
        returnvalue += CalcFTGDiagonalTerm(meas_y - (ul_corner_y + y_size),
            meas_z - (ul_corner_z + z_size), meas_x - (ul_corner_x + x_size));

        return returnvalue * Grav_const;
      }

    /*! Calculate the geometric term of the xy-element of the second order derivative tensor, we use the same
     * function with permutated parameters to calculate the other off-diagonal elements.
     * The description of the parameters is for the Uxy calculation.
     * @param meas_x x-coordinate of the measurement in m
     * @param meas_y y-coordinate of the measurement in m
     * @param meas_z z-coordinate of the measurement in m
     * @param ul_corner_x x-coordinate of the upper left front corner of the density element in m
     * @param ul_corner_y y-coordinate of the upper left front corner of the density element in m
     * @param ul_corner_z z-coordinate of the upper left front corner of the density element in m
     * @param x_size size of the prism in x-direction in m
     * @param y_size size of the prism in y-direction in m
     * @param z_size size of the prism in z-direction in m
     * @return The Uxx element of the FTG tensor
     */
    double CalcUxyTerm(const double meas_x, const double meas_y,
        const double meas_z, const double ul_corner_x,
        const double ul_corner_y, const double ul_corner_z,
        const double x_size, const double y_size, const double z_size)
      {
        //we need 8 terms, one for each corner of the prism, they have varying signs
        double returnvalue = -CalcFTGOffDiagonalTerm(meas_z - ul_corner_z,
            meas_y - ul_corner_y, meas_z - ul_corner_z, meas_x - ul_corner_x);
        // one corner shifted -> positive sign
        returnvalue += CalcFTGOffDiagonalTerm(meas_z - ul_corner_z, meas_y
            - ul_corner_y, meas_z - ul_corner_z, meas_x
            - (ul_corner_x + x_size));
        returnvalue += CalcFTGOffDiagonalTerm(meas_z - ul_corner_z, meas_y
            - (ul_corner_y + y_size), meas_z - ul_corner_z, meas_x
            - ul_corner_x);
        returnvalue += CalcFTGOffDiagonalTerm(meas_z - (ul_corner_z + z_size),
            meas_y - ul_corner_y, meas_z - (ul_corner_z + z_size), meas_x
                - ul_corner_x);
        // two corners shifted -> negative sign
        returnvalue -= CalcFTGOffDiagonalTerm(meas_z - ul_corner_z, meas_y
            - (ul_corner_y + y_size), meas_z - ul_corner_z, meas_x
            - (ul_corner_x + x_size));
        returnvalue -= CalcFTGOffDiagonalTerm(meas_z - (ul_corner_z + z_size),
            meas_y - ul_corner_y, meas_z - (ul_corner_z + z_size), meas_x
                - (ul_corner_x + x_size));
        returnvalue -= CalcFTGOffDiagonalTerm(meas_z - (ul_corner_z + z_size),
            meas_y - (ul_corner_y + y_size), meas_z - (ul_corner_z + z_size),
            meas_x - ul_corner_x);
        //three corners shifted -> positive sign
        returnvalue += CalcFTGOffDiagonalTerm(meas_z - (ul_corner_z + z_size),
            meas_y - (ul_corner_y + y_size), meas_z - (ul_corner_z + z_size),
            meas_x - (ul_corner_x + x_size));

        return -returnvalue * Grav_const;
      }

    /*! Calculate FTG tensor for a rectangular prism
     * @param meas_x x-coordinate of the measurement in m
     * @param meas_y y-coordinate of the measurement in m
     * @param meas_z z-coordinate of the measurement in m
     * @param ul_corner_x x-coordinate of the upper left front corner of the density element in m
     * @param ul_corner_y y-coordinate of the upper left front corner of the density element in m
     * @param ul_corner_z z-coordinate of the upper left front corner of the density element in m
     * @param x_size size of the prism in x-direction in m
     * @param y_size size of the prism in y-direction in m
     * @param z_size size of the prism in z-direction in m
     * @return The tensor response of the prism at the measurement site
     */
    rmat CalcTensorBoxTerm(const double meas_x,
        const double meas_y, const double meas_z, const double ul_corner_x,
        const double ul_corner_y, const double ul_corner_z,
        const double x_size, const double y_size, const double z_size)
      {
        rmat returnvalue(3, 3);
        //Calculate the first diagonal element
        returnvalue(0, 0) = CalcUxxTerm(meas_x, meas_y, meas_z, ul_corner_x,
            ul_corner_y, ul_corner_z, x_size, y_size, z_size);
        // The other two diagonal elements can be calculated from permutations of the above equation
        returnvalue(1, 1) = CalcUxxTerm(meas_y, meas_x, meas_z, ul_corner_y,
            ul_corner_x, ul_corner_z, y_size, x_size, z_size);
        returnvalue(2, 2) = CalcUxxTerm(meas_z, meas_x, meas_y, ul_corner_z,
            ul_corner_x, ul_corner_y, z_size, x_size, y_size);
        //Calculate the first off-diagonal element
        returnvalue(0, 1) = CalcUxyTerm(meas_x, meas_y, meas_z, ul_corner_x,
            ul_corner_y, ul_corner_z, x_size, y_size, z_size);
        //the gravimetry matrix is symmetric
        returnvalue(1, 0) = returnvalue(0, 1);
        //The other off-diagonal can be calculated from permutations of the above equation
        returnvalue(0, 2) = CalcUxyTerm(meas_x, meas_z, meas_y, ul_corner_x,
            ul_corner_z, ul_corner_y, x_size, z_size, y_size);
        returnvalue(2, 0) = returnvalue(0, 2);
        returnvalue(1, 2) = CalcUxyTerm(meas_z, meas_y, meas_x, ul_corner_z,
            ul_corner_y, ul_corner_x, z_size, y_size, x_size);
        returnvalue(2, 1) = returnvalue(1, 2);

        return returnvalue;
      }

    /*! Calculate the vertical accelerational effect of a semi-infinite sheet with constant density
     * @param hor_dist Horizontal distance of the sheet from the measurement site in m
     * @param ver_dist Vertical distance of the sheet from the measurement site in m
     * @param thick Thickness of the sheet in m
     * @param density Density of the sheet in g/cm^3
     * @return Gravitational acceleration in m/s^2
     */
    double CalcGravSemiInfSheet(const double hor_dist, const double ver_dist,
        const double thick, const double density)
      {
        return (2.0 * Grav_const * density * thick) * ((M_PI / 2.0) - atan2(
            hor_dist, ver_dist));
      }

  }
