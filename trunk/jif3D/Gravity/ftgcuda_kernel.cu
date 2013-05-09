#ifndef _FTGCUDA_KERNEL_H_
#define _FTGCUDA_KERNEL_H_

#include "common.cu"

__device__ double CalcFTGDiagonalTerm(const double a, const double b,
    const double c)
  {
    return atan2(a * b, c * sqrt(a * a + b * b + c * c));
  }

__device__ double CalcFTGOffDiagonalTerm(const double value, const double x,
    const double y, const double z)
  {
    return log(value + sqrt(x * x + y * y + z * z));
  }

__device__ double CalcUxxTerm(const double meas_x, const double meas_y,
    const double meas_z, const double ul_corner_x, const double ul_corner_y,
    const double ul_corner_z, const double x_size, const double y_size,
    const double z_size)
  {
    //we need 8 terms, one for each corner of the prism, they have varying signs
    double returnvalue = -CalcFTGDiagonalTerm(meas_y - ul_corner_y, meas_z
        - ul_corner_z, meas_x - ul_corner_x);
    // one corner shifted -> positive sign
    returnvalue += CalcFTGDiagonalTerm(meas_y - ul_corner_y, meas_z
        - ul_corner_z, meas_x - (ul_corner_x + x_size));
    returnvalue += CalcFTGDiagonalTerm(meas_y - (ul_corner_y + y_size), meas_z
        - ul_corner_z, meas_x - ul_corner_x);
    returnvalue += CalcFTGDiagonalTerm(meas_y - ul_corner_y, meas_z
        - (ul_corner_z + z_size), meas_x - ul_corner_x);
    // two corners shifted -> negative sign
    returnvalue -= CalcFTGDiagonalTerm(meas_y - (ul_corner_y + y_size), meas_z
        - ul_corner_z, meas_x - (ul_corner_x + x_size));
    returnvalue -= CalcFTGDiagonalTerm(meas_y - ul_corner_y, meas_z
        - (ul_corner_z + z_size), meas_x - (ul_corner_x + x_size));
    returnvalue -= CalcFTGDiagonalTerm(meas_y - (ul_corner_y + y_size), meas_z
        - (ul_corner_z + z_size), meas_x - ul_corner_x);
    //three corners shifted -> positive sign
    returnvalue += CalcFTGDiagonalTerm(meas_y - (ul_corner_y + y_size), meas_z
        - (ul_corner_z + z_size), meas_x - (ul_corner_x + x_size));

    return returnvalue * Grav_const;
  }
//! Calculate the Uxy geometric term for a single rectangular prism, is also used for the other off-diagonal elements with permuted arguments
__device__ double CalcUxyTerm(const double meas_x, const double meas_y,
    const double meas_z, const double ul_corner_x, const double ul_corner_y,
    const double ul_corner_z, const double x_size, const double y_size,
    const double z_size)
  {
    //we need 8 terms, one for each corner of the prism, they have varying signs
    double returnvalue = -CalcFTGOffDiagonalTerm(meas_z - ul_corner_z, meas_y
        - ul_corner_y, meas_z - ul_corner_z, meas_x - ul_corner_x);
    // one corner shifted -> positive sign
    returnvalue += CalcFTGOffDiagonalTerm(meas_z - ul_corner_z, meas_y
        - ul_corner_y, meas_z - ul_corner_z, meas_x - (ul_corner_x + x_size));
    returnvalue += CalcFTGOffDiagonalTerm(meas_z - ul_corner_z, meas_y
        - (ul_corner_y + y_size), meas_z - ul_corner_z, meas_x - ul_corner_x);
    returnvalue += CalcFTGOffDiagonalTerm(meas_z - (ul_corner_z + z_size),
        meas_y - ul_corner_y, meas_z - (ul_corner_z + z_size), meas_x
            - ul_corner_x);
    // two corners shifted -> negative sign
    returnvalue -= CalcFTGOffDiagonalTerm(meas_z - ul_corner_z, meas_y
        - (ul_corner_y + y_size), meas_z - ul_corner_z, meas_x - (ul_corner_x
        + x_size));
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

__global__ void CalcUxxMeas(const double x_meas, const double y_meas,
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
        //currvalue contains the geometric term, i.e. the sensitivity

        returnvalue[offset] = CalcUxxTerm(x_meas, y_meas, z_meas,
            XCoord[xindex], YCoord[yindex], ZCoord[zindex], XSizes[xindex],
            YSizes[yindex], ZSizes[zindex]);

      }
  }

__global__ void CalcUyyMeas(const double x_meas, const double y_meas,
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
        //currvalue contains the geometric term, i.e. the sensitivity

        returnvalue[offset] = CalcUxxTerm(y_meas, x_meas, z_meas,
            YCoord[yindex], XCoord[xindex], ZCoord[zindex], YSizes[yindex],
            XSizes[xindex], ZSizes[zindex]);

      }
  }

__global__ void CalcUzzMeas(const double x_meas, const double y_meas,
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
        //currvalue contains the geometric term, i.e. the sensitivity

        returnvalue[offset] = CalcUxxTerm(z_meas, x_meas, y_meas,
            ZCoord[zindex], XCoord[xindex], YCoord[yindex], ZSizes[zindex],
            XSizes[xindex], YSizes[yindex]);

      }
  }


__global__ void CalcUxyMeas(const double x_meas, const double y_meas,
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
        //currvalue contains the geometric term, i.e. the sensitivity

        returnvalue[offset] = CalcUxyTerm(x_meas, y_meas, z_meas,
            XCoord[xindex], YCoord[yindex], ZCoord[zindex], XSizes[xindex],
            YSizes[yindex], ZSizes[zindex]);

      }
  }

__global__ void CalcUxzMeas(const double x_meas, const double y_meas,
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
        //currvalue contains the geometric term, i.e. the sensitivity

        returnvalue[offset] = CalcUxyTerm(x_meas, z_meas, y_meas,
            XCoord[xindex], ZCoord[zindex], YCoord[yindex], XSizes[xindex],
            ZSizes[zindex], YSizes[yindex]);

      }
  }

__global__ void CalcUyzMeas(const double x_meas, const double y_meas,
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
        //currvalue contains the geometric term, i.e. the sensitivity

        returnvalue[offset] = CalcUxyTerm(z_meas, y_meas, x_meas,
            ZCoord[zindex], YCoord[yindex], XCoord[xindex], ZSizes[zindex],
            YSizes[yindex], XSizes[xindex]);

      }
  }


#endif
