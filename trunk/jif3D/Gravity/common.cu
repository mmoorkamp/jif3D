#ifndef _COMMON_H_
#define _COMMON_H_


#define Grav_const 6.67428e-11

__device__ void OffsetToIndex(const unsigned int offset, const int ny,
    const int nz, int &xi, int &yi, int &zi)
  {
    zi = offset % nz;
    xi = (offset - zi) / nz;
    yi = xi % ny;
    xi = (xi - yi) / ny;
  }


#endif
