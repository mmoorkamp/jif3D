//============================================================================
// Name        : ThreeDGravityModel.cpp
// Author      : Max Moorkamp
// Version     :
// Copyright   : 2008, MM
//============================================================================

#include "ThreeDGravityModel.h"
#include "../Global/NumUtil.h"
#include "ReadWriteGravityData.h"
#include <algorithm>
#include <boost/bind.hpp>
#include <cassert>
#include <netcdfcpp.h>
#include <fstream>
#include <iomanip>
#include <cmath>

namespace jiba
  {
    static const std::string DensityName = "density";
    static const std::string DensityUnit = "g/cm3";

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
     * @param ul_corner_x x-coordinate of the upper left front corner of the density element
     * @param ul_corner_y y-coordinate of the upper left front corner of the density element
     * @param ul_corner_z z-coordinate of the upper left front corner of the density element
     * @param x_size size of the prism in x-direction in m
     * @param y_size size of the prism in y-direction in m
     * @param z_size size of the prism in z-direction in m
     * @return The tensor response of the prism at the measurement site
     */
    GravimetryMatrix CalcTensorBoxTerm(const double meas_x,
        const double meas_y, const double meas_z, const double ul_corner_x,
        const double ul_corner_y, const double ul_corner_z,
        const double x_size, const double y_size, const double z_size)
      {
        GravimetryMatrix returnvalue(3, 3);
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
    // Calculate a single term of the equation for a  diagonal  FTG element
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
                    //we store the current value for possible sensitivity calculations
                    currvalue = CalcGravBoxTerm(x_meas, y_meas, z_meas,
                        GetXCoordinates()[i], GetYCoordinates()[j],
                        GetZCoordinates()[k], GetXCellSizes()[i],
                        GetYCellSizes()[j], GetZCellSizes()[k]);
                    returnvalue += currvalue * GetDensities()[i][j][k];
                    // if we want to store the sensitivity matrix
                    if (StoreScalarSensitivities)
                      {
                        //we calculate the sensitivity by dividing through density
                        //therefore density cannot be zero
                        ScalarSensitivities(meas_index, i * (ysize * zsize) + j
                            * zsize + k) = currvalue;
                      }
                  }
              }
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
                    currvalue = CalcTensorBoxTerm(x_meas, y_meas, z_meas,
                        GetXCoordinates()[i], GetYCoordinates()[j],
                        GetZCoordinates()[k], GetXCellSizes()[i],
                        GetYCellSizes()[j], GetZCellSizes()[k]);
                    returnvalue += currvalue * GetDensities()[i][j][k];
                    //if we want to store sensitivities for tensor measurements
                    if (StoreTensorSensitivities)
                      {
                        boost::numeric::ublas::matrix_column<rmat> column(
                            TensorSensitivities, i * (ysize * zsize) + j
                                * zsize + k); // extract the right column of the sensitivity matrix
                        std::copy(currvalue.data().begin(),
                            currvalue.data().end(), // assign the elements to the right part of the column
                            column.begin() + meas_index * 9);
                      }
                  }
              }
          }

        return returnvalue;
      }

    /*!  Calculate the contribution of a layered background to a scalar gravity measurement.
     * @param xmeas The x-coordinate of the measurement point in m
     * @param ymeas The y-coordinate of the measurement point in m
     * @param zmeas The z-coordinate of the measurement point in m
     * @param xwidth The total width of the discretized model area in x-direction in m
     * @param ywidth The total width of the discretized model area in y-direction in m
     * @param zwidth The total width of the discretized model area in z-direction in m
     * @param meas_index The index of the measurement
     * @return The gravitational acceleration in m/s^2 due to the background
     */
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
            currvalue = CalcInfSheetTerm(bg_thicknesses[j]);
            currbottom = currtop + bg_thicknesses[j];
            // and then subtract the value for the modelling domain, as this is already calculated in the discretized routine
            // if the background layer complete coincides with the discretized area
            if (currtop < zwidth && (currbottom <= zwidth))

              {
                currvalue -= CalcGravBoxTerm(xmeas, ymeas, zmeas, 0.0, 0.0,
                    currtop, xwidth, ywidth, bg_thicknesses[j]);
              }
            //if some of the background coincides and some is below
            if (currtop < zwidth && currbottom > zwidth)

              {
                currvalue -= CalcGravBoxTerm(xmeas, ymeas, zmeas, 0.0, 0.0,
                    currtop, xwidth, ywidth, (zwidth - currtop));
              }
            //we also store the sensitivities for the background
            if (StoreScalarSensitivities)
              {
                ScalarSensitivities(meas_index, modelsize + j) = currvalue;
              }
            result += currvalue * bg_densities[j];
            currtop += bg_thicknesses[j];
          }
        return result;
      }

    /*!  Calculate the contribution of a layered background to a tensor gravity measurement.
     * @param xmeas The x-coordinate of the measurement point in m
     * @param ymeas The y-coordinate of the measurement point in m
     * @param zmeas The z-coordinate of the measurement point in m
     * @param xwidth The total width of the discretized model area in x-direction in m
     * @param ywidth The total width of the discretized model area in y-direction in m
     * @param zwidth The total width of the discretized model area in z-direction in m
     * @param meas_index The index of the measurement
     * @return The gravitational tensor due to the background
     */
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
                currvalue = -CalcTensorBoxTerm(x_meas, y_meas, z_meas, 0.0,
                    0.0, currtop, xwidth, ywidth, bg_thicknesses[j]);

              }
            if (currtop < zwidth && currbottom > zwidth) //if some of the background coincides and some is below

              {
                currvalue = -CalcTensorBoxTerm(x_meas, y_meas, z_meas, 0.0,
                    0.0, currtop, xwidth, ywidth, (zwidth - currtop));
              }
            if (StoreTensorSensitivities)
              {
                boost::numeric::ublas::matrix_column<rmat> column(
                    TensorSensitivities, modelsize + j); // extract the right column of the sensitivity matrix
                std::copy(currvalue.data().begin(), currvalue.data().end(), // assign the elements to the right part of the column
                    column.begin() + meas_index * 9);
              }
            result += currvalue * bg_densities[j];
            currtop += bg_thicknesses[j];
          }
        return result;
      }

    ThreeDGravityModel::tScalarMeasVec ThreeDGravityModel::CalcGravity()
      {
        //get the amount of cells in each direction
        const size_t xsize = GetData().shape()[0];
        const size_t ysize = GetData().shape()[1];
        const size_t zsize = GetData().shape()[2];
        //do some sanity checks
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
        ScalarResults.assign(nmeas, 0.0);

        // if we have already stored the ScalarSensitivity we can do it very fast
        if (HaveCalculatedScalarSensitivities)
          {
            //we just construct a vector of densities
            rvec DensityVector(xsize * ysize * zsize + bg_densities.size());
            copy(GetDensities().origin(), GetDensities().origin()
                + GetDensities().num_elements(), DensityVector.begin());

            //including the background
            copy(bg_densities.begin(), bg_densities.end(),
                DensityVector.begin() + xsize * ysize * zsize);
            //and do a Matrix Vector multiplication
            rvec MeasVec(prec_prod(ScalarSensitivities, DensityVector));
            copy(MeasVec.begin(), MeasVec.end(), ScalarResults.begin());
            return ScalarResults;
          }
        // we only get here if we didn't store the sensitivities previously
        //check if we want to store sensitivities and allocate memory
        if (StoreScalarSensitivities && !HaveCalculatedScalarSensitivities)
          {
            ScalarSensitivities.resize(MeasPosX.size(), xsize * ysize * zsize
                + bg_densities.size());
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
            ScalarResults[i] = CalcScalarMeas(MeasPosX[i], MeasPosY[i],
                MeasPosZ[i], i);
            ScalarResults[i] += CalcBackground(MeasPosX[i], MeasPosY[i],
                MeasPosZ[i], modelxwidth, modelywidth, modelzwidth, i);

          }
        // if we store the sensitivity matrix
        if (StoreScalarSensitivities)
          {
            //we remember that we have it so we can accelerate the next calculations
            HaveCalculatedScalarSensitivities = true;
          }
        return ScalarResults;
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
        const size_t nmeas = MeasPosX.size();
        assert(nmeas == MeasPosY.size());
        // make sure we have coordinates for all sites
        assert(nmeas == MeasPosZ.size());
        // calculate the size of the modelling domain for the background adjustment
        const double modelxwidth = std::accumulate(GetXCellSizes().begin(),
            GetXCellSizes().end(), 0.0);
        const double modelywidth = std::accumulate(GetYCellSizes().begin(),
            GetYCellSizes().end(), 0.0);
        const double modelzwidth = std::accumulate(GetZCellSizes().begin(),
            GetZCellSizes().end(), 0.0);

        TensorResults.assign(nmeas, rmat(3, 3));
        // if we have already stored the ScalarSensitivity we can do it very fast
        if (HaveCalculatedTensorSensitivities)
          {
            //create a vector of discretized densities and background densities
            rvec DensityVector(xsize * ysize * zsize + bg_densities.size());
            copy(GetDensities().origin(), GetDensities().origin()
                + GetDensities().num_elements(), DensityVector.begin());
            copy(bg_densities.begin(), bg_densities.end(),
                DensityVector.begin() + xsize * ysize * zsize);
            //do a vector matrix multiplication
            rvec MeasVec(prec_prod(TensorSensitivities, DensityVector));
            //and copy the resulting vector in the right format for the return value
            for (size_t i = 0; i < nmeas; ++i)
              {
                copy(MeasVec.begin() + i * 9, MeasVec.begin() + (i + 1) * 9,
                    TensorResults.at(i).data().begin());
              }
            return TensorResults;
          }
        //we only get here if we didn't use the stored sensitivities
        if (StoreTensorSensitivities && !HaveCalculatedTensorSensitivities)
          {
            TensorSensitivities.resize(nmeas * 9, xsize * ysize * zsize
                + bg_densities.size()); // we have 9 tensor elements for each measurement points
          }
        // for all measurement points add the responses of the discretized part and the 1D background
        for (size_t i = 0; i < nmeas; ++i)
          {
            TensorResults[i] = CalcTensorMeas(MeasPosX[i], MeasPosY[i],
                MeasPosZ[i], i);
            //adjust for the effects of finite extents of the grid
            TensorResults[i] += AdjustTensorBackground(MeasPosX[i],
                MeasPosY[i], MeasPosZ[i], modelxwidth, modelywidth,
                modelzwidth, i);

          }
        if (StoreTensorSensitivities)
          {
            //remember that we have already calculated these values
            HaveCalculatedTensorSensitivities = true;
          }
        return TensorResults;

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
        WriteDataToNetCDF(DataFile, DensityName, DensityUnit);
      }

    void ThreeDGravityModel::ReadNetCDF(const std::string filename)
      {
        NcFile DataFile(filename.c_str(), NcFile::ReadOnly);
        ReadDataFromNetCDF(DataFile, DensityName, DensityUnit);
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

    void ThreeDGravityModel::SaveScalarMeasurements(const std::string filename)
      {
        SaveScalarGravityMeasurements(filename, ScalarResults, MeasPosX,
            MeasPosY, MeasPosZ);
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
        PlotMeasAscii(filename, ScalarResults);
      }

    void ThreeDGravityModel::PlotTensorMeasurements(
        const std::string filename_root)
      {
        //at the moment we only write out the diagonal elements
        tScalarMeasVec Uxx, Uyy, Uzz;
        const size_t nmeas = TensorResults.size();
        for (size_t i = 0; i < nmeas; ++i)
          {
            Uxx.push_back(TensorResults.at(i)(0, 0));
            Uyy.push_back(TensorResults.at(i)(1, 1));
            Uzz.push_back(TensorResults.at(i)(2, 2));
          }
        PlotMeasAscii(filename_root + ".uxx.plot", Uxx);
        PlotMeasAscii(filename_root + ".uyy.plot", Uyy);
        PlotMeasAscii(filename_root + ".uzz.plot", Uzz);
      }

    void ThreeDGravityModel::ReadMeasPosNetCDF(const std::string filename)
      {
        jiba::ReadMeasPosNetCDF(filename, MeasPosX, MeasPosY, MeasPosZ);
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
