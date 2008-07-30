//============================================================================
// Name        : ThreeDGravityModel.cpp
// Author      : Max Moorkamp
// Version     :
// Copyright   : 2008, MM
//============================================================================

#include "ThreeDGravityModel.h"
#include "../Global/NumUtil.h"
#include <iostream>
#include <algorithm>
#include <boost/bind.hpp>
#include <cassert>
#include <netcdfcpp.h>
#include <fstream>
#include <iomanip>
#include <cmath>

namespace jiba
  {
    static const std::string GravDataName = "density";
    static const std::string GravDataUnit = "g/cm^3";

    //! Calculate one term for the gravitational potential of a box, we use the nomenclature of eq. 4-6 in Li and Chouteau
    /*! Calculate one term for the gravitational potential of a box, we use the nomenclature of eq. 4-6 in Li and Chouteau.
     * The parameters x,y and z are the distances to the corners of the box, we will call this functions with different
     * permutations of real world x,y and z coordinates, so only in some cases x corresponds to the x-axis.
     */
    double CalcGravTerm(const double x, const double y, const double z)
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

    //! Calculate the gravitational potential of a rectangular prism at a point meas_{x,y,z}
    /*! Calculate the gravitational potential of a rectangular prism at a point meas_{x,y,z}
     *  The calculation is based on equation 4 in Li and Chouteau, Surveys in Geophysics, 19, 339-368, 1998
     * @param meas_x x-coordinate of the measurement in m
     * @param meas_y y-coordinate of the measurement in m
     * @param meas_z z-coordinate of the measurement in m
     * @param ul_corner_x x-coordinate of the upper left front corner of the density element
     * @param ul_corner_y y-coordinate of the upper left front corner of the density element
     * @param ul_corner_z z-coordinate of the upper left front corner of the density element
     * @param x_size size of the prism in x-direction in m
     * @param y_size size of the prism in y-direction in m
     * @param z_size size of the prism in z-direction in m
     * @param density density of the prism in g/cm^3
     * @return Gravitational acceleration in m/s^2
     */
    double CalcGravBox(const double meas_x, const double meas_y,
        const double meas_z, const double ul_corner_x,
        const double ul_corner_y, const double ul_corner_z,
        const double x_size, const double y_size, const double z_size,
        const double density)
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
        return -returnvalue * Grav_const * density;
      }
    //! Calculate the xx-element of the second order derivative tensor
    /*! Calculate the xx-element of the second order derivative tensor, we use the same
     * function with permutated parameters to calculate the other two diagonal elements.
     * The description of the parameters is for the Uxx calculation.
     * @param meas_x x-coordinate of the measurement in m
     * @param meas_y y-coordinate of the measurement in m
     * @param meas_z z-coordinate of the measurement in m
     * @param ul_corner_x x-coordinate of the upper left front corner of the density element
     * @param ul_corner_y y-coordinate of the upper left front corner of the density element
     * @param ul_corner_z z-coordinate of the upper left front corner of the density element
     * @param x_size size of the prism in x-direction in m
     * @param y_size size of the prism in y-direction in m
     * @param z_size size of the prism in z-direction in m
     * @param density density of the prism in g/cm^3
     * @return The Uxx element of the FTG tensor
     */
    double CalcUxx(const double meas_x, const double meas_y,
        const double meas_z, const double ul_corner_x,
        const double ul_corner_y, const double ul_corner_z,
        const double x_size, const double y_size, const double z_size,
        const double density)
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

        return returnvalue * Grav_const * density;
      }
    //! Calculate the xy-element of the second order derivative tensor
    /*! Calculate the xy-element of the second order derivative tensor, we use the same
     * function with permutated parameters to calculate the other off-diagonal elements.
     * The description of the parameters is for the Uxy calculation.
     * @param meas_x x-coordinate of the measurement in m
     * @param meas_y y-coordinate of the measurement in m
     * @param meas_z z-coordinate of the measurement in m
     * @param ul_corner_x x-coordinate of the upper left front corner of the density element
     * @param ul_corner_y y-coordinate of the upper left front corner of the density element
     * @param ul_corner_z z-coordinate of the upper left front corner of the density element
     * @param x_size size of the prism in x-direction in m
     * @param y_size size of the prism in y-direction in m
     * @param z_size size of the prism in z-direction in m
     * @param density density of the prism in g/cm^3
     * @return The Uxx element of the FTG tensor
     */
    double CalcUxy(const double meas_x, const double meas_y,
        const double meas_z, const double ul_corner_x,
        const double ul_corner_y, const double ul_corner_z,
        const double x_size, const double y_size, const double z_size,
        const double density)
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

        return -returnvalue * Grav_const * density;
      }
    //! Calculate the second order derivative tensor for a box shaped element
    /*! Calculate FTG tensor for a rectangular prism
     * @param meas_x x-coordinate of the measurement in m
     * @param meas_y y-coordinate of the measurement in m
     * @param meas_z z-coordinate of the measurement in m
     * @param ul_corner_x x-coordinate of the upper left front corner of the density element
     * @param ul_corner_y y-coordinate of the upper left front corner of the density element
     * @param ul_corner_z z-coordinate of the upper left front corner of the density element
     * @param x_size size of the prism in x-direction in m
     * @param y_size size of the prism in y-direction in m
     * @param z_size size of the prism in z-direction in m
     * @param density density of the prism in g/cm^3
     * @return The tensor response of the prism at the measurement site
     */
    GravimetryMatrix CalcTensorBox(const double meas_x, const double meas_y,
        const double meas_z, const double ul_corner_x,
        const double ul_corner_y, const double ul_corner_z,
        const double x_size, const double y_size, const double z_size,
        const double density)
      {
        GravimetryMatrix returnvalue(3, 3);
        //Calculate the first diagonal element
        returnvalue(0, 0) = CalcUxx(meas_x, meas_y, meas_z, ul_corner_x,
            ul_corner_y, ul_corner_z, x_size, y_size, z_size, density);
        // The other two diagonal elements can be calculated from permutations of the above equation
        returnvalue(1, 1) = CalcUxx(meas_y, meas_x, meas_z, ul_corner_y,
            ul_corner_x, ul_corner_z, y_size, x_size, z_size, density);
        returnvalue(2, 2) = CalcUxx(meas_z, meas_x, meas_y, ul_corner_z,
            ul_corner_x, ul_corner_y, z_size, x_size, y_size, density);
        //Calculate the first off-diagonal element
        returnvalue(0, 1) = CalcUxy(meas_x, meas_y, meas_z, ul_corner_x,
            ul_corner_y, ul_corner_z, x_size, y_size, z_size, density);
        //the gravimetry matrix is symmetric
        returnvalue(1, 0) = returnvalue(0, 1);
        //The other off-diagonal can be calculated from permutations of the above equation
        returnvalue(0, 2) = CalcUxy(meas_x, meas_z, meas_y, ul_corner_x,
            ul_corner_z, ul_corner_y, x_size, z_size, y_size, density);
        returnvalue(2, 0) = returnvalue(0, 2);
        returnvalue(1, 2) = CalcUxy(meas_z, meas_y, meas_x, ul_corner_z,
            ul_corner_y, ul_corner_x, z_size, y_size, x_size, density);
        returnvalue(2, 1) = returnvalue(1, 2);

        return returnvalue;
      }
    //! Calculate the effect of a semi-infinite sheet
    /*! Calculate the effect of a semi-infinite sheet with constant density
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
    //! Calculate a single term of the equation for a  diagonal  FTG element
    double CalcFTGDiagonalTerm(const double a, const double b, const double c)
      {
        return atan2(a * b, c * sqrt(a * a + b * b + c * c));
        // this is an alternative formula to check the correctness
        //        const double r = sqrt(a*a+b*b+c*c);
        //
        //        double returnvalue = atan2(a * b, c *r);
        //        returnvalue += a*b*c*r/(c*c*r*r + a*a*b*b);
        //        returnvalue += pow(c,3) * b * a/(c*c*pow(r,3)+b*b*a*a*r);
        //        returnvalue += a*c/(b*r+r*r) + b*c/(a*r+r*r);
        //        return returnvalue;
      }
    //! Calculate a single term of the equation for an off-diagonal FTG element
    double CalcFTGOffDiagonalTerm(const double value, const double x,
        const double y, const double z)
      {
        return log(value + std::max(sqrt(x * x + y * y + z * z),
            std::numeric_limits<double>::epsilon()));
      }

    //! Calculate the gravitational effect of the 3D model at a single measurement site
    /*! Calculate the gravitational effect of the 3D model at a single measurement site. The way we calculate
     * the sensitivity matrix at the moment, the model cannot contain densities of 0 if we
     * store the sensitivity matrix
     * @param x_meas x-coordinate of the measurement site in m
     * @param y_meas y-coordinate of the measurement site in m
     * @param z_meas z-coordinate of the measurement site in m
     * @param meas_index index of the measurement site among all measurements, this is only for storing sensitivities
     * @return The gravitational acceleration in m/s^2 due to the model at this site
     */
    double ThreeDGravityModel::CalcScalarMeas(const double x_meas,
        const double y_meas, const double z_meas, const size_t meas_index)
      {
        //get the dimensions of the model
        const size_t xsize = GetData().shape()[0];
        const size_t ysize = GetData().shape()[1];
        const size_t zsize = GetData().shape()[2];

        double returnvalue = 0.0;
        double currvalue = 0.0;
        //sum up the contributions of all prisms
        for (size_t i = 0; i < xsize; ++i)
          {
            for (size_t j = 0; j < ysize; ++j)
              {
                for (size_t k = 0; k < zsize; ++k)
                  {
                    currvalue = CalcGravBox(x_meas, y_meas, z_meas,
                        GetXCoordinates()[i], GetYCoordinates()[j],
                        GetZCoordinates()[k], GetXCellSizes()[i],
                        GetYCellSizes()[j], GetZCellSizes()[k],
                        GetDensities()[i][j][k]);
                    returnvalue += currvalue;
                    // if we want to store the sensitivity matrix
                    if (StoreScalarSensitivities)
                      {
                        //we calculate the sensitivity by dividing through density
                        //therefore density cannot be zero
                        ScalarSensitivities(meas_index, i * (ysize * zsize) + j
                            * zsize + k) = currvalue / GetDensities()[i][j][k];
                      }
                  }
              }
          }
        // if we store the sensitivity matrix
        if (StoreScalarSensitivities)
          {
            //we remember that we have it so we can accelerate the next calculations
            HaveCalculatedScalarSensitivities = true;
          }
        return returnvalue;
      }

    GravimetryMatrix ThreeDGravityModel::CalcTensorMeas(const double x_meas,
        const double y_meas, const double z_meas, const size_t meas_index)
      {
        const size_t xsize = GetData().shape()[0];
        const size_t ysize = GetData().shape()[1];
        const size_t zsize = GetData().shape()[2];

        GravimetryMatrix returnvalue(3, 3);
        returnvalue *= 0.0;
        GravimetryMatrix currvalue(3, 3);
        //sum up the contributions of all prisms
        for (size_t i = 0; i < xsize; ++i)
          {
            for (size_t j = 0; j < ysize; ++j)
              {
                for (size_t k = 0; k < zsize; ++k)
                  {
                    currvalue = CalcTensorBox(x_meas, y_meas, z_meas,
                        GetXCoordinates()[i], GetYCoordinates()[j],
                        GetZCoordinates()[k], GetXCellSizes()[i],
                        GetYCellSizes()[j], GetZCellSizes()[k],
                        GetDensities()[i][j][k]);
                    returnvalue += currvalue;
                    if (StoreTensorSensitivities)
                      {
                        boost::numeric::ublas::matrix_column<rmat> column(
                            TensorSensitivities, i * (ysize * zsize) + j
                                * zsize + k); // extract the right column of the sensitivity matrix
                        std::transform(currvalue.data().begin(),
                            currvalue.data().end(), // assign the elements to the right part of the column
                            column.begin() + meas_index * 9, boost::bind(
                                std::divides<double>(), _1,
                                GetDensities()[i][j][k]));
                      }
                  }
              }
          }
        if (StoreScalarSensitivities)
          {
            HaveCalculatedScalarSensitivities = true;
          }
        return returnvalue;
      }

    double ThreeDGravityModel::CalcBackground(const double xmeas,
        const double ymeas, const double zmeas, const double xwidth,
        const double ywidth, const double zwidth, const size_t meas_index)
      {
        //make sure we have thicknesses and densities for all layers
        assert(bg_densities.size() == bg_thicknesses.size());
        const size_t nbglayers = bg_densities.size();
        double result = 0.0;
        double currtop = 0.0;
        double currvalue = 0.0;
        double currbottom = 0.0;
        const size_t modelsize = GetData().shape()[0] * GetData().shape()[1]
            * GetData().shape()[2];
        // for all layers of the background
        for (size_t j = 0; j < nbglayers; ++j)
          {
            // first assume an infinite sheet for the current layer
            currvalue = CalcInfSheet(bg_thicknesses[j], bg_densities[j]);
            currbottom = currtop + bg_thicknesses[j];
            // and then subtract the value for the modelling domain, as this is already calculated in the discretized routine
            // if the background layer complete coincides with the discretized area
            if (currtop < zwidth && (currbottom <= zwidth))

              {
                currvalue -= CalcGravBox(xmeas, ymeas, zmeas, 0.0, 0.0, currtop,
                        xwidth, ywidth, bg_thicknesses[j], bg_densities[j]);
              }
            //if some of the background coincides and some is below
            if (currtop < zwidth && currbottom > zwidth)

              {
                currvalue -= CalcGravBox(xmeas, ymeas, zmeas, 0.0, 0.0,
                    currtop, xwidth, ywidth, (zwidth - currtop),
                    bg_densities[j]);
              }
            if (StoreScalarSensitivities)
              {
                ScalarSensitivities(meas_index, modelsize + j) = currvalue
                    / bg_densities[j];
              }
            result += currvalue;
            currtop += bg_thicknesses[j];
          }
        return result;
      }

    GravimetryMatrix ThreeDGravityModel::AdjustTensorBackground(
        const double x_meas, const double y_meas, const double z_meas,
        const double xwidth, const double ywidth, const double zwidth,
        const size_t meas_index)
      {
        //make sure we have thicknesses and densities for all layers
        assert(bg_densities.size() == bg_thicknesses.size());
        const size_t nbglayers = bg_densities.size();
        GravimetryMatrix result(3, 3);
        GravimetryMatrix currvalue(3, 3);
        result *= 0.0;
        double currtop = 0.0;
        double currbottom = 0.0;
        const size_t modelsize = GetData().shape()[0] * GetData().shape()[1]
            * GetData().shape()[2];
        // for all layers of the background
        for (size_t j = 0; j < nbglayers; ++j)
          {
            currbottom = currtop + bg_thicknesses[j];
            if (currtop < zwidth && (currbottom <= zwidth)) // if the background layer complete coincides with the discretized area

              {
                // We only have to substract the effect of the gridding box, the effect of an inifite sheet is zero
                currvalue
                    = -CalcTensorBox(x_meas, y_meas, z_meas, 0.0, 0.0, currtop,
                        xwidth, ywidth, bg_thicknesses[j], bg_densities[j]);
              }
            if (currtop < zwidth && currbottom > zwidth) //if some of the background coincides and some is below

              {
                currvalue = -CalcTensorBox(x_meas, y_meas, z_meas, 0.0, 0.0,
                    currtop, xwidth, ywidth, (zwidth - currtop),
                    bg_densities[j]);
              }
            if (StoreTensorSensitivities)
              {
                boost::numeric::ublas::matrix_column<rmat> column(
                    TensorSensitivities, modelsize + j); // extract the right column of the sensitivity matrix
                std::transform(currvalue.data().begin(),
                    currvalue.data().end(), // assign the elements to the right part of the column
                    column.begin() + meas_index * 9, boost::bind(
                        std::divides<double>(), _1, bg_densities[j]));
              }
            result += currvalue;
            currtop += bg_thicknesses[j];
          }
        return result;
      }

    ThreeDGravityModel::tScalarMeasVec ThreeDGravityModel::CalcGravity()
      {
        const size_t xsize = GetData().shape()[0];
        const size_t ysize = GetData().shape()[1];
        const size_t zsize = GetData().shape()[2];

        assert(xsize == GetXCoordinates().shape()[0]);
        assert(ysize == GetYCoordinates().shape()[0]);
        assert(zsize == GetZCoordinates().shape()[0]);
        assert(xsize == GetXCellSizes().shape()[0]);
        assert(ysize == GetYCellSizes().shape()[0]);
        assert(zsize == GetZCellSizes().shape()[0]);

        // make sure we have coordinates for all sites
        const size_t nmeas = MeasPosX.size();
        assert(nmeas == MeasPosY.size());
        assert(nmeas == MeasPosZ.size());
        tScalarMeasVec results(nmeas);

        // if we have already stored the ScalarSensitivity we can do it very fast
        if (HaveCalculatedScalarSensitivities)
          {
            rvec DensityVector(xsize * ysize * zsize + bg_densities.size());
            for (size_t i = 0; i < xsize; ++i)
              {
                for (size_t j = 0; j < ysize; ++j)
                  {
                    for (size_t k = 0; k < zsize; ++k)
                      {
                        DensityVector(i * (ysize * zsize) + j * zsize + k)
                            = GetDensities()[i][j][k];
                      }
                  }
              }
            copy(bg_densities.begin(), bg_densities.end(),
                DensityVector.begin() + xsize * ysize * zsize);

            rvec MeasVec(prec_prod(ScalarSensitivities, DensityVector));
            copy(MeasVec.begin(), MeasVec.end(), results.begin());
            return results;
          }
        // we only get here if we didn't store the sensitivities previously
        //check if store sensitivities and allocate memory
        if (StoreScalarSensitivities && !HaveCalculatedScalarSensitivities)
          {
            ScalarSensitivities.resize(MeasPosX.size(), xsize * ysize * zsize
                + bg_densities.size());
          }
        if (StoreTensorSensitivities && !HaveCalculatedTensorSensitivities)
          {
            TensorSensitivities.resize(nmeas * 9, xsize * ysize * zsize); // we have 9 tensor elements for each measurement points
          }
        // calculate the size of the modelling domain for the background adjustment
        const double modelxwidth = std::accumulate(GetXCellSizes().begin(),
            GetXCellSizes().end(), 0.0);
        const double modelywidth = std::accumulate(GetYCellSizes().begin(),
            GetYCellSizes().end(), 0.0);
        const double modelzwidth = std::accumulate(GetZCellSizes().begin(),
            GetZCellSizes().end(), 0.0);

        // for all measurement points add the respones of the discretized part and the 1D background
        for (size_t i = 0; i < nmeas; ++i)
          {
            results[i] = CalcScalarMeas(MeasPosX[i], MeasPosY[i], MeasPosZ[i],
                i);
            results[i] += CalcBackground(MeasPosX[i], MeasPosY[i], MeasPosZ[i],
                modelxwidth, modelywidth, modelzwidth, i);

          }
        return results;
      }

    ThreeDGravityModel::tTensorMeasVec ThreeDGravityModel::CalcTensorGravity()
      {
        const size_t xsize = GetData().shape()[0];
        const size_t ysize = GetData().shape()[1];
        const size_t zsize = GetData().shape()[2];

        assert(xsize == GetXCoordinates().shape()[0]);
        assert(ysize == GetYCoordinates().shape()[0]);
        assert(zsize == GetZCoordinates().shape()[0]);
        assert(xsize == GetXCellSizes().shape()[0]);
        assert(ysize == GetYCellSizes().shape()[0]);
        assert(zsize == GetZCellSizes().shape()[0]);

        assert(MeasPosX.size() == MeasPosY.size());
        // make sure we have coordinates for all sites
        assert(MeasPosX.size() == MeasPosZ.size());
        // calculate the size of the modelling domain for the background adjustment
        const double modelxwidth = std::accumulate(GetXCellSizes().begin(),
            GetXCellSizes().end(), 0.0);
        const double modelywidth = std::accumulate(GetYCellSizes().begin(),
            GetYCellSizes().end(), 0.0);
        const double modelzwidth = std::accumulate(GetZCellSizes().begin(),
            GetZCellSizes().end(), 0.0);

        const size_t nmeas = MeasPosX.size();
        tTensorMeasVec results(nmeas, rmat(3, 3));
        // if we have already stored the ScalarSensitivity we can do it very fast
        if (HaveCalculatedTensorSensitivities)
          {
            rvec DensityVector(xsize * ysize * zsize + bg_densities.size());
            for (size_t i = 0; i < xsize; ++i)
              {
                for (size_t j = 0; j < ysize; ++j)
                  {
                    for (size_t k = 0; k < zsize; ++k)
                      {
                        DensityVector(i * (ysize * zsize) + j * zsize + k)
                            = GetDensities()[i][j][k];
                      }
                  }
              }
            copy(bg_densities.begin(), bg_densities.end(),
                DensityVector.begin() + xsize * ysize * zsize);

            rvec MeasVec(prec_prod(ScalarSensitivities, DensityVector));
            for (size_t i = 0; i < nmeas; ++i)
              {
                copy(MeasVec.begin() + i * 9, MeasVec.begin() + (i + 1) * 9,
                    results.at(i).data().begin());
              }
            return results;
          }
        // for all measurement points add the respones of the discretized part and the 1D background
        for (size_t i = 0; i < nmeas; ++i)
          {
            results[i] = CalcTensorMeas(MeasPosX[i], MeasPosY[i], MeasPosZ[i],
                i);
            //adjust for the effects of finite extents of the grid
            results[i] += AdjustTensorBackground(MeasPosX[i], MeasPosY[i],
                MeasPosZ[i], modelxwidth, modelywidth, modelzwidth, i);

          }
        return results;

      }

    ThreeDGravityModel::ThreeDGravityModel(const bool storescalar,
        const bool storetensor) :
      StoreScalarSensitivities(storescalar), StoreTensorSensitivities(
          storetensor), HaveCalculatedScalarSensitivities(false),
          HaveCalculatedTensorSensitivities(false)
      {
      }

    ThreeDGravityModel::~ThreeDGravityModel()
      {
      }

    void ThreeDGravityModel::WriteNetCDF(const std::string filename) const
      {
        NcFile DataFile(filename.c_str(), NcFile::Replace);
        WriteDataToNetCDF(DataFile, GravDataName, GravDataUnit);
      }

    void ThreeDGravityModel::ReadNetCDF(const std::string filename)
      {
        NcFile DataFile(filename.c_str(), NcFile::ReadOnly);
        ReadDataFromNetCDF(DataFile, GravDataName, GravDataUnit);
      }

    NcDim *ThreeDGravityModel::WriteDimensionToNetCDF(NcFile &NetCDFFile,
        const std::string &SizeName, const tMeasPosVec &Position) const
      {

        // Add a dimension and a variable with the same name to the netcdf file
        NcDim *SizeDim = NetCDFFile.add_dim(SizeName.c_str(), Position.size());
        NcVar *SizeVar =
            NetCDFFile.add_var(SizeName.c_str(), ncDouble, SizeDim);
        //All length is measured in meters
        SizeVar->add_att("units", "m");
        //We also store the name
        SizeVar->add_att("long_name", (SizeName + " coordinate").c_str());
        // we store the coordinates of the cells in the netcdf file

        //Write the values
        SizeVar->put(&Position[0], Position.size());
        // We return the NcDim object, because we need it to write the model data
        return SizeDim;
      }

    void ThreeDGravityModel::ReadDimensionFromNetCDF(NcFile &NetCDFFile,
        const std::string &DimName, tMeasPosVec &Position)
      {

        //create a netcdf dimension with the chosen name
        NcDim *Dim = NetCDFFile.get_dim(DimName.c_str());
        //determine the size of that dimension
        const size_t nvalues = Dim->size();

        //allocate memory in the class variable
        Position.assign(nvalues, 0.0);
        // create netcdf variable with the same name as the dimension
        NcVar *SizeVar = NetCDFFile.get_var(DimName.c_str());
        //read coordinate values from netcdf file
        SizeVar->get(&Position[0], nvalues);
      }

    void ThreeDGravityModel::SaveScalarMeasurements(const std::string filename)
      {
        tScalarMeasVec Data = CalcGravity();

        assert(Data.size() == MeasPosX.size());
        assert(Data.size() == MeasPosY.size());
        assert(Data.size() == MeasPosZ.size());
        const size_t nmeas = MeasPosX.size();
        NcFile DataFile(filename.c_str(), NcFile::Replace);
        std::vector<int> StationNumber;
        std::generate_n(back_inserter(StationNumber), nmeas, IntSequence(0));
        NcDim *StatNumDim = DataFile.add_dim("StationNumber", nmeas);
        NcVar *StatNumVar =
            DataFile.add_var("StationNumber", ncInt, StatNumDim);
        StatNumVar->put(&StationNumber[0], nmeas);
        NcVar *XPosVar = DataFile.add_var("xpos", ncDouble, StatNumDim);
        XPosVar->add_att("units", "m");
        XPosVar->put(&MeasPosX[0], nmeas);
        NcVar *YPosVar = DataFile.add_var("ypos", ncDouble, StatNumDim);
        YPosVar->add_att("units", "m");
        YPosVar->put(&MeasPosY[0], nmeas);
        NcVar *ZPosVar = DataFile.add_var("zpos", ncDouble, StatNumDim);
        ZPosVar->add_att("units", "m");
        ZPosVar->put(&MeasPosZ[0], nmeas);
        //Write the measurements
        NcVar *DataVar = DataFile.add_var("Scalar_gravity", ncDouble,
            StatNumDim);
        DataVar->add_att("units", "m/s^2");
        DataVar->add_att("_FillValue", -1.0);

        DataVar->put(&Data[0], StatNumDim->size());
      }

    void ThreeDGravityModel::PlotMeasAscii(const std::string &filename,
        tScalarMeasVec &Data) const
      {
        std::ofstream outfile(filename.c_str());

        const size_t nmeas = Data.size();
        assert(nmeas == MeasPosX.size());
        assert(nmeas == MeasPosY.size());
        assert(nmeas == MeasPosZ.size());
        for (size_t i = 0; i < nmeas; ++i)
          {
            outfile << std::setw(15) << std::setprecision(5) << MeasPosX.at(i)
                << " ";
            outfile << std::setw(15) << std::setprecision(5) << MeasPosY.at(i)
                << " ";
            outfile << std::setw(15) << std::setprecision(5) << MeasPosZ.at(i)
                << " ";
            outfile << std::setw(15) << std::setprecision(5) << Data.at(i);
            outfile << std::endl;
          }
      }

    void ThreeDGravityModel::PlotScalarMeasurements(const std::string filename)
      {
        tScalarMeasVec Data(CalcGravity());

        PlotSingleMeasAscii(filename, Data);

        //    const size_t nmeas = MeasPosX.size();
        //        assert(Data.size() == MeasPosX.size());
        //        assert(Data.size() == MeasPosY.size());
        //
        //
        //        NcFile DataFile(filename.c_str(), NcFile::Replace);
        //        // Write the size information in x,y, and z-direction
        //        NcDim *XSizeDim = WriteDimensionToNetCDF(DataFile, "x", MeasPosX);
        //        NcDim *YSizeDim = WriteDimensionToNetCDF(DataFile, "y", MeasPosY);
        //
        //        NcVar *DataVar = DataFile.add_var("Scalar_gravity", ncDouble,
        //                    XSizeDim, YSizeDim);
        //        DataVar->add_att("units", "m/s^2");
        //        DataVar->add_att("_FillValue",-1.0);
        //        //Write the measurements
        //        tScalarMeasVec WrittenData(nmeas*nmeas,-1.0);
        //            for (size_t i = 0; i < nmeas; ++i)
        //              {
        //              WrittenData.at(i + nmeas*i ) = Data.at(i);
        //              }
        //            DataVar->put(&WrittenData[0], XSizeDim->size(), YSizeDim->size());
      }

    void ThreeDGravityModel::PlotTensorMeasurements(
        const std::string filename_root)
      {
        tTensorMeasVec Data(CalcTensorGravity());

        tScalarMeasVec Uxx, Uyy, Uzz;
        const size_t nmeas = Data.size();
        for (size_t i = 0; i < nmeas; ++i)
          {
            Uxx.push_back(Data.at(i)(0, 0));
            Uyy.push_back(Data.at(i)(1, 1));
            Uzz.push_back(Data.at(i)(2, 2));
          }
        PlotSingleMeasAscii(filename_root + ".uxx.plot", Uxx);
        PlotSingleMeasAscii(filename_root + ".uyy.plot", Uyy);
        PlotSingleMeasAscii(filename_root + ".uzz.plot", Uzz);
      }

    void ThreeDGravityModel::ReadMeasPosNetCDF(const std::string filename)
      {
        NcFile DataFile(filename.c_str(), NcFile::ReadOnly);
        ReadDimensionFromNetCDF(DataFile, "Measx", MeasPosX);
        ReadDimensionFromNetCDF(DataFile, "Measy", MeasPosY);
        ReadDimensionFromNetCDF(DataFile, "Measz", MeasPosZ);
        assert(MeasPosX.size() == MeasPosY.size());
        assert(MeasPosX.size() == MeasPosZ.size());
      }

    void ThreeDGravityModel::ReadMeasPosAscii(const std::string filename)
      {
        std::ifstream infile(filename.c_str());
        double posx, posy, posz;
        while (infile.good())
          {
            infile >> posx >> posy >> posz;
            if (infile.good())
              {
                MeasPosX.push_back(posx);
                MeasPosY.push_back(posy);
                MeasPosZ.push_back(posz);
              }
          }
        assert(MeasPosX.size() == MeasPosY.size());
        assert(MeasPosX.size() == MeasPosZ.size());
      }
  }
