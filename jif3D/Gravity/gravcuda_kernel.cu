#ifndef _GRAVCUDA_KERNEL_H_
#define _GRAVCUDA_KERNEL_H_

#include "common.cu"
//The calculation of the effect of a rectilinear prism involves 8 terms
//that all have the same form, but with different permutations of the arguments
//this function implements these elementary terms
__device__ double CalcGravTerm(const double x, const double y, const double z)
  {
    const double r = sqrt(x * x + y * y + z * z);
    double rvalue = x * log(y + r) + y * log(x + r);
    // atan2 takes care of small denominators
    rvalue += z * atan2(z * r, x * y);
    return rvalue;
  }
//calculate the geometric term of the effect of a rectilinear grid element
//on a scalar gravity measurement station
__device__ double CalcGravBoxTerm(const double meas_x, const double meas_y,
    const double meas_z, const double ul_corner_x, const double ul_corner_y,
    const double ul_corner_z, const double x_size, const double y_size,
    const double z_size)
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
    // we multiply the geometric term by the gravitational constant
    //the multiplication by density is performed by the host
    //this gives us more flexibility for the calculation of the gradient
    return -returnvalue * Grav_const;
  }
//the implementation of the CUDA kernel, calculates the geometric terms for a 3D grid
//of rectilinear mesh elements on a single gravity measurement station
__global__ void CalcScalarMeas(const double x_meas, const double y_meas,
    const double z_meas, const double *XCoord, const double *YCoord,
    const double *ZCoord, const double *XSizes, const double *YSizes,
    const double *ZSizes, const int nx, const int ny, const int nz,
    double *returnvalue)
  {
    const unsigned int offset = blockIdx.x * blockDim.x + threadIdx.x;
    int xindex, yindex, zindex;
    if (offset < nx * ny * nz)
      {
        OffsetToIndex(offset, ny, nz, xindex, yindex, zindex);
        //we store the current value for possible sensitivity calculations
        //returnvalue contains the geometric term, i.e. the sensitivity
        returnvalue[offset] = CalcGravBoxTerm(x_meas, y_meas, z_meas,
            XCoord[xindex], YCoord[yindex], ZCoord[zindex], XSizes[xindex],
            YSizes[yindex], ZSizes[zindex]);

      }
  }

#endif // #ifndef _GRAVCUDA_KERNEL_H_
