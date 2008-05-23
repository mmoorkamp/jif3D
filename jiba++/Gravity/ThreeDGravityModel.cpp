#include "ThreeDGravityModel.h"
#include "../Global/NumUtil.h"
#include <iostream>
#include <algorithm>
#include <boost/bind.hpp>

namespace jiba
  {
    static const std::string GravDataName = "density";
    static const std::string GravDataUnit = "g/cm^3";

    //! Calculate one term for the gravitational potential of a box, we use the nomenclature of eq. 4-6 in Li and Chouteau
    double CalcGravTerm(const double x, const double y, const double z)
      {
        //assert(fabs(x*y)> std::numeric_limits<double>::epsilon());
        const double r = sqrt(x*x+y*y+z*z);
        double rvalue = x * log(y +r) + y * log(x+r);
        if (fabs(x*y)> std::numeric_limits<double>::epsilon())
          {
            rvalue += z*atan2(z*r, x*y);
          }
        else
          {
            if (sign(x*y) == sign(z*r))
              rvalue += z* M_PI/2.0;
            else
              rvalue -= z* M_PI/2.0;
          }

        return rvalue;
      }

    //! Calculate the gravitational potential of a rectangular prism at a point meas_{x,y,z}
    /*! The calculation is based on equation 4 in Li and Chouteau, Surveys in Geophysics, 19, 339-368, 1998
     */
    double CalcGravBox(const double meas_x, const double meas_y,
        const double meas_z, const double ul_corner_x,
        const double ul_corner_y, const double ul_corner_z,
        const double x_size, const double y_size, const double z_size,
        const double density)
      {
        //we need 8 terms, one for each corner of the prism, they have varying signs
        double returnvalue = -CalcGravTerm(meas_x-ul_corner_x, meas_y
            -ul_corner_y, meas_z-ul_corner_z);
        // one corner shifted -> positive sign
        returnvalue += CalcGravTerm(meas_x-(ul_corner_x+x_size), meas_y
            -ul_corner_y, meas_z-ul_corner_z);
        returnvalue += CalcGravTerm(meas_x-ul_corner_x, meas_y -(ul_corner_y
            +y_size), meas_z-ul_corner_z);
        returnvalue += CalcGravTerm(meas_x-ul_corner_x, meas_y -ul_corner_y,
            meas_z-(ul_corner_z+z_size));
        // two corners shifted -> negative sign
        returnvalue -= CalcGravTerm(meas_x-(ul_corner_x+x_size), meas_y
            -(ul_corner_y+y_size), meas_z-ul_corner_z);
        returnvalue -= CalcGravTerm(meas_x-(ul_corner_x+x_size), meas_y
            -ul_corner_y, meas_z-(ul_corner_z+z_size));
        returnvalue -= CalcGravTerm(meas_x-ul_corner_x, meas_y -(ul_corner_y
            +y_size), meas_z-(ul_corner_z+z_size));
        //three corners shifted -> positive sign
        returnvalue += CalcGravTerm(meas_x-(ul_corner_x+x_size), meas_y
            -(ul_corner_y+y_size), meas_z-(ul_corner_z+z_size));

        return -returnvalue * Grav_const * density;
      }

    double CalcUxx(const double meas_x, const double meas_y,
        const double meas_z, const double ul_corner_x,
        const double ul_corner_y, const double ul_corner_z,
        const double x_size, const double y_size, const double z_size,
        const double density)
      {
        //we need 8 terms, one for each corner of the prism, they have varying signs
        double returnvalue = -CalcFTGDiagonalTerm(meas_y -ul_corner_y, meas_z
            -ul_corner_z, meas_x-ul_corner_x);
        // one corner shifted -> positive sign
        returnvalue += CalcFTGDiagonalTerm(meas_y -ul_corner_y, meas_z
            -ul_corner_z, meas_x-(ul_corner_x+x_size));
        returnvalue += CalcFTGDiagonalTerm(meas_y -(ul_corner_y +y_size),
            meas_z -ul_corner_z, meas_x-ul_corner_x);
        returnvalue += CalcFTGDiagonalTerm(meas_y -ul_corner_y, meas_z
            -(ul_corner_z +z_size), meas_x-ul_corner_x);
        // two corners shifted -> negative sign
        returnvalue -= CalcFTGDiagonalTerm(meas_y -(ul_corner_y+y_size), meas_z
            -ul_corner_z, meas_x-(ul_corner_x+x_size));
        returnvalue -= CalcFTGDiagonalTerm(meas_y -ul_corner_y, meas_z
            -(ul_corner_z +z_size), meas_x-(ul_corner_x+x_size));
        returnvalue -= CalcFTGDiagonalTerm(meas_y -(ul_corner_y +y_size),
            meas_z -(ul_corner_z+z_size), meas_x-ul_corner_x);
        //three corners shifted -> positive sign
        returnvalue += CalcFTGDiagonalTerm(meas_y -(ul_corner_y+y_size), meas_z
            -(ul_corner_z+z_size), meas_x-(ul_corner_x+x_size));

        return returnvalue * Grav_const * density;
      }

    double CalcUxy(const double meas_x, const double meas_y,
        const double meas_z, const double ul_corner_x,
        const double ul_corner_y, const double ul_corner_z,
        const double x_size, const double y_size, const double z_size,
        const double density)
      {
        //we need 8 terms, one for each corner of the prism, they have varying signs
        double returnvalue = -CalcFTGOffDiagonalTerm(meas_z -ul_corner_z,
            meas_y -ul_corner_y, meas_z -ul_corner_z, meas_x-ul_corner_x);
        // one corner shifted -> positive sign
        returnvalue += CalcFTGOffDiagonalTerm(meas_z -ul_corner_z, meas_y
            -ul_corner_y, meas_z -ul_corner_z, meas_x-(ul_corner_x+x_size));
        returnvalue += CalcFTGOffDiagonalTerm(meas_z -ul_corner_z, meas_y
            -(ul_corner_y +y_size), meas_z -ul_corner_z, meas_x-ul_corner_x);
        returnvalue += CalcFTGOffDiagonalTerm(meas_z -(ul_corner_z +z_size),
            meas_y -ul_corner_y, meas_z -(ul_corner_z +z_size), meas_x
                -ul_corner_x);
        // two corners shifted -> negative sign
        returnvalue -= CalcFTGOffDiagonalTerm(meas_z -ul_corner_z, meas_y
            -(ul_corner_y+y_size), meas_z -ul_corner_z, meas_x-(ul_corner_x
            +x_size));
        returnvalue -= CalcFTGOffDiagonalTerm(meas_z -(ul_corner_z +z_size),
            meas_y -ul_corner_y, meas_z -(ul_corner_z +z_size), meas_x
                -(ul_corner_x+x_size));
        returnvalue -= CalcFTGOffDiagonalTerm(meas_z -(ul_corner_z+z_size),
            meas_y -(ul_corner_y +y_size), meas_z -(ul_corner_z+z_size), meas_x
                -ul_corner_x);
        //three corners shifted -> positive sign
        returnvalue += CalcFTGOffDiagonalTerm(meas_z -(ul_corner_z+z_size),
            meas_y -(ul_corner_y+y_size), meas_z -(ul_corner_z+z_size), meas_x
                -(ul_corner_x+x_size));

        return -returnvalue * Grav_const * density;
      }

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

    double CalcGravSemiInfSheet(const double hor_dist, const double ver_dist,
        const double thick, const double density)
      {
        return (2.0*Grav_const*density*thick)*((M_PI/2.0) - atan2(hor_dist,
            ver_dist));
      }

    double CalcFTGDiagonalTerm(const double a, const double b, const double c)
      {
        return -atan2(a * b, c *sqrt(a*a+b*b+c*c));
      }

    double CalcFTGOffDiagonalTerm(const double value, const double x,
        const double y, const double z)
      {
        return log(value + sqrt(x*x+y*y+z*z));
      }

    double ThreeDGravityModel::CalcScalarMeas(const double x_meas,
        const double y_meas, const double z_meas, const size_t meas_index)
      {
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
                        GetXCoordinates()[i], GetYCoordinates()[j], GetZCoordinates()[k], GetXCellSizes()[i], GetYCellSizes()[j], GetZCellSizes()[k], GetDensities()[i][j][k]);
                    returnvalue += currvalue;
                    if (StoreScalarSensitivities)
                      {
                        ScalarSensitivities(meas_index, i*(ysize*zsize)+j*zsize
                            +k) = currvalue/GetDensities()[i][j][k];
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
                        GetXCoordinates()[i], GetYCoordinates()[j], GetZCoordinates()[k], GetXCellSizes()[i], GetYCellSizes()[j], GetZCellSizes()[k], GetDensities()[i][j][k]);
                    returnvalue += currvalue;
                    if (StoreTensorSensitivities)
                      {
                        boost::numeric::ublas::matrix_column<rmat> column(
                            TensorSensitivities, i*(ysize*zsize)+j*zsize +k); // extract the right column of the sensitivity matrix
                        std::transform(currvalue.data().begin(), currvalue.data().end(), // assign the elements to the right part of the column 
                            column.begin()+meas_index*9, boost::bind(
                                std::divides<double>(),_1,GetDensities()[i][j][k]));
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
                const size_t modelsize = GetData().shape()[0] * GetData().shape()[1] * GetData().shape()[2];
        // for all layers of the background
        for (size_t j = 0; j < nbglayers; ++j)
          {
            // first assume an infinite sheet for the current layer
            currvalue = CalcInfSheet(bg_thicknesses[j], bg_densities[j]);
            currbottom = currtop + bg_thicknesses[j];
            // and then subtract the value for the modelling domain, as this is already calculated in the discretized routine
            if (currtop < zwidth && (currbottom <= zwidth)) // if the background layer complete coincides with the discretized area

              {
                currvalue
                -= CalcGravBox(xmeas, ymeas, zmeas, 0.0, 0.0, currtop,
                xwidth, ywidth, bg_thicknesses[j], bg_densities[j]);
              }
            if (currtop < zwidth && currbottom> zwidth) //if some of the background coincides and some is below

              {
                currvalue
                -= CalcGravBox(xmeas, ymeas, zmeas, 0.0, 0.0, currtop,
                xwidth, ywidth, (zwidth - currtop), bg_densities[j]);
              }
            if (StoreScalarSensitivities)
              {
                ScalarSensitivities(meas_index, modelsize+j) = currvalue
                /bg_densities[j];
              }
            result += currvalue;
            currtop += bg_thicknesses[j];
          }
        return result;
      }

    GravimetryMatrix ThreeDGravityModel::AdjustTensorBackground(
    const double x_meas, const double y_meas, const double z_meas,
    const double xwidth, const double ywidth, const double zwidth, const size_t meas_index)
      {
        //make sure we have thicknesses and densities for all layers
        assert(bg_densities.size() == bg_thicknesses.size());
        const size_t nbglayers = bg_densities.size();
        GravimetryMatrix result(3, 3);
        GravimetryMatrix currvalue(3, 3);
        result *= 0.0;
        double currtop = 0.0;
        double currbottom = 0.0;
        const size_t modelsize = GetData().shape()[0] * GetData().shape()[1] * GetData().shape()[2];
        // for all layers of the background
        for (size_t j = 0; j < nbglayers; ++j)
          {
            currbottom = currtop + bg_thicknesses[j];
            if (currtop < zwidth && (currbottom <= zwidth)) // if the background layer complete coincides with the discretized area

              {
                // We only have to substract the effect of the gridding box, the effect of an inifite sheet is zero
                currvalue = -CalcTensorBox(x_meas, y_meas, z_meas, 0.0, 0.0, currtop,
                xwidth, ywidth, bg_thicknesses[j], bg_densities[j]);
              }
            if (currtop < zwidth && currbottom> zwidth) //if some of the background coincides and some is below

              {
                currvalue = -CalcTensorBox(x_meas, y_meas, z_meas, 0.0, 0.0, currtop,
                xwidth, ywidth, (zwidth - currtop), bg_densities[j]);
              }
            if (StoreTensorSensitivities)
              {
                boost::numeric::ublas::matrix_column<rmat> column(
                TensorSensitivities, modelsize + j); // extract the right column of the sensitivity matrix
                std::transform(currvalue.data().begin(), currvalue.data().end(), // assign the elements to the right part of the column 
                column.begin()+meas_index*9, boost::bind(
                std::divides<double>(),_1,bg_densities[j]));
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
        rvec DensityVector(xsize*ysize*zsize+bg_densities.size());
        for (size_t i = 0; i < xsize; ++i)
          {
            for (size_t j = 0; j < ysize; ++j)
              {
                for (size_t k = 0; k < zsize; ++k)
                  {
                    DensityVector(i*(ysize*zsize)+j*zsize+k) = GetDensities()[i][j][k];
                  }
              }
          }
        copy(bg_densities.begin(), bg_densities.end(), DensityVector.begin()
            +xsize*ysize*zsize);

        rvec MeasVec(prec_prod(ScalarSensitivities, DensityVector));
        copy(MeasVec.begin(), MeasVec.end(), results.begin());
        return results;
      }
    // we only get here if we didn't store the sensitivities previously
    //check if store sensitivities and allocate memory
    if (StoreScalarSensitivities && !HaveCalculatedScalarSensitivities)
      {
        ScalarSensitivities.resize(MeasPosX.size(), xsize*ysize*zsize
            +bg_densities.size());
      }
    if (StoreTensorSensitivities && !HaveCalculatedTensorSensitivities)
      {
        TensorSensitivities.resize(nmeas * 9, xsize*ysize*zsize); // we have 9 tensor elements for each measurement points
      }
    // calculate the size of the modelling domain for the background adjustment
    const double modelxwidth = std::accumulate(GetXCellSizes().begin(), GetXCellSizes().end(), 0.0);
    const double modelywidth = std::accumulate(GetYCellSizes().begin(), GetYCellSizes().end(), 0.0);
    const double modelzwidth = std::accumulate(GetZCellSizes().begin(), GetZCellSizes().end(), 0.0);

    // for all measurement points add the respones of the discretized part and the 1D background
    for (size_t i = 0; i < nmeas; ++i)
      {
        results[i] = CalcScalarMeas(MeasPosX[i], MeasPosY[i], MeasPosZ[i], i);
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
    const double modelxwidth = std::accumulate(GetXCellSizes().begin(), GetXCellSizes().end(), 0.0);
    const double modelywidth = std::accumulate(GetYCellSizes().begin(), GetYCellSizes().end(), 0.0);
    const double modelzwidth = std::accumulate(GetZCellSizes().begin(), GetZCellSizes().end(), 0.0);

    const size_t nmeas = MeasPosX.size();
    tTensorMeasVec results(nmeas, rmat(3, 3));
    // if we have already stored the ScalarSensitivity we can do it very fast
    if (HaveCalculatedTensorSensitivities)
      {
        rvec DensityVector(xsize*ysize*zsize+bg_densities.size());
        for (size_t i = 0; i < xsize; ++i)
          {
            for (size_t j = 0; j < ysize; ++j)
              {
                for (size_t k = 0; k < zsize; ++k)
                  {
                    DensityVector(i*(ysize*zsize)+j*zsize+k) = GetDensities()[i][j][k];
                  }
              }
          }
        copy(bg_densities.begin(), bg_densities.end(), DensityVector.begin()
            +xsize*ysize*zsize);

        rvec MeasVec(prec_prod(ScalarSensitivities, DensityVector));
        for (size_t i = 0; i < nmeas; ++i)
          {
            copy(MeasVec.begin()+i*9, MeasVec.begin()+(i+1)*9, results.at(i).data().begin());
          }
        return results;
      }
    // for all measurement points add the respones of the discretized part and the 1D background
    for (size_t i = 0; i < nmeas; ++i)
      {
        results[i] = CalcTensorMeas(MeasPosX[i], MeasPosY[i], MeasPosZ[i], i);
        //adjust for the effects of finite extents of the grid
        results[i] += AdjustTensorBackground(MeasPosX[i], MeasPosY[i],
            MeasPosZ[i], modelxwidth, modelywidth, modelzwidth, i);

      }
    return results;

  }

ThreeDGravityModel::ThreeDGravityModel(const bool storescalar,
    const bool storetensor) :
  StoreScalarSensitivities(storescalar), StoreTensorSensitivities(storetensor),
      HaveCalculatedScalarSensitivities(false),
      HaveCalculatedTensorSensitivities(false)
  {
  }

ThreeDGravityModel::~ThreeDGravityModel()
  {
  }

void ThreeDGravityModel::WriteNetCDF(const std::string filename)
  {
    NcFile DataFile(filename.c_str(), NcFile::Replace);
    WriteDataToNetCDF(DataFile, GravDataName, GravDataUnit);
  }

void ThreeDGravityModel::ReadNetCDF(const std::string filename)
  {
    NcFile DataFile(filename.c_str(), NcFile::ReadOnly);
    ReadDataFromNetCDF(DataFile, GravDataName, GravDataUnit);
  }

}
