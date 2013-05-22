//============================================================================
// Name        : Wavelet.cpp
// Author      : Sep 26, 2008
// Version     :
// Copyright   : 2008, mmoorkamp
//============================================================================

#include "Wavelet.h"
#include "FatalException.h"
#include "NumUtil.h"
#include <numeric>

/*! \file Wavelet.cpp
 * Simple implementations of the discrete wavelet transform and its inverse based on Daubechies-4 wavelets.
 * The basic algorithm is described in Numerical Recipes in C, pp. 595-598
 */

namespace jif3D
  {
    //calculate the coefficients for the wavelet filter
    static const double c0 = (1.0 + std::sqrt(3.0)) / (4.0 * std::sqrt(2.0));
    static const double c1 = (3.0 + std::sqrt(3.0)) / (4.0 * std::sqrt(2.0));
    static const double c2 = (3.0 - std::sqrt(3.0)) / (4.0 * std::sqrt(2.0));
    static const double c3 = (1.0 - std::sqrt(3.0)) / (4.0 * std::sqrt(2.0));

    /*! Apply a 4 coefficient Daubechies wavelet filter to the input vector.
     * @param Invec The vector with the input data, will contain the result of the transform. Size must be a power of 2
     * @param maxindex The maximum index to which to apply the filter, must be a power of 2
     */
    void Daub4(jif3D::rvec &Invec, const size_t maxindex)
      {
        //we perform the operations on a temporary copy
        jif3D::rvec Temp(Invec);
        const size_t length = maxindex;
        const size_t offset = length / 2;
        //the current value depends on the next three values, so we have to adjust the end of the loop
        //also note the different stride for j and i
        for (size_t i = 0, j = 0; j < length - 3; j += 2, ++i)
          {
            Temp(i) = c0 * Invec(j) + c1 * Invec(j + 1) + c2 * Invec(j + 2)
                + c3 * Invec(j + 3);
            Temp(i + offset) = c3 * Invec(j) - c2 * Invec(j + 1) + c1 * Invec(j
                + 2) - c0 * Invec(j + 3);
          }
        //we have to handle the last values as a special case
        Temp(length / 2 - 1) = c0 * Invec(length - 2) + c1 * Invec(length - 1)
            + c2 * Invec(0) + c3 * Invec(1);
        Temp(length - 1) = c3 * Invec(length - 2) - c2 * Invec(length - 1) + c1
            * Invec(0) - c0 * Invec(1);
        //copy the result back to Invec, we only copy up to maxindex
        std::copy(Temp.begin(), Temp.begin() + maxindex, Invec.begin());
      }


    /*! Apply a 4 coefficient inverse Daubechies wavelet filter to the input vector.
     * @param Invec The vector with the input data, will contain the result of the transform. Size must be a power of 2
     * @param maxindex The maximum index to which to apply the filter, must be a power of 2 and less than the size of Invec
     */
    void InvDaub4(jif3D::rvec &Invec, const size_t maxindex)
      {
        //this is the same strategy as the forward
        //only the equations and the order are slightly different
        jif3D::rvec Temp(Invec);
        const size_t length = maxindex;
        const size_t offset = length / 2;
        Temp(0) = c2 * Invec(length / 2 - 1) + c1 * Invec(length - 1) + c0
            * Invec(0) + c3 * Invec(length / 2);
        Temp(1) = c3 * Invec(length / 2 - 1) - c0 * Invec(length - 1) + c1
            * Invec(0) - c2 * Invec(length / 2);
        for (size_t i = 0, j = 2; i < offset - 1; ++i, j += 2)
          {
            Temp(j) = c2 * Invec(i) + c1 * Invec(i + offset) + c0
                * Invec(i + 1) + c3 * Invec(i + offset + 1);
            Temp(j + 1) = c3 * Invec(i) - c0 * Invec(i + offset) + c1 * Invec(i
                + 1) - c2 * Invec(i + offset + 1);
          }
        std::copy(Temp.begin(), Temp.begin() + maxindex, Invec.begin());
      }

    /*! The classes ForwardWaveletDriver and InverseWaveletDriver capture the
     * code that is different for the forward and inverse transforms in one
     * and multidimensional transformations. They are basically the one-dimensional transforms.
     */
    class ForwardWaveletDriver
      {
    public:
      /*! The function operator takes two arguments
       * @param WorkVector The vector with the input data, will contain the transformed data
       * @param DimSize The maximum length the transform will be performed on
       */
      void operator()(rvec &WorkVector, size_t DimSize)
        {
          for (size_t n = DimSize; n >= 4; n >>= 1)
            {
              Daub4(WorkVector, n);
            }
        }
      };

    class InverseWaveletDriver
      {
    public:
      /*! The function operator takes two arguments
       * @param WorkVector The vector with the input data, will contain the transformed data
       * @param DimSize The maximum length the transform will be performed on
       */
      void operator()(rvec &WorkVector, size_t DimSize)
        {
          for (size_t n = 4; n <= DimSize; n <<= 1)
            {
              InvDaub4(WorkVector, n);
            }
        }
      };

    /*! Perform a forward wavelet transform on the input vector
     * @param Invec The input vector, will contain the result, size must be power of 2
     */
    void WaveletTransform(jif3D::rvec &Invec)
      {
        ForwardWaveletDriver()(Invec, Invec.size());
      }

    /*! Perform an inverse wavelet transform on the input vector
     * @param Invec The input vector, will contain the result, size must be power of 2
     */
    void InvWaveletTransform(jif3D::rvec &Invec)
      {
        InverseWaveletDriver()(Invec, Invec.size());
      }

    /*! The structure of the multidimensional wavelet transform is identical
     * for forward and inverse. We therefore put the difference in a separate
     * function object and pass this as a parameter. For effectivity and because
     * we know at compile time whether its inverse or forward we make it a template.
     * @param InArray The array with the input values, contains the output afterwards
     * @param DimSizes The size of each dimension, InArray must contain the product of DimSizes elements
     * @param ndim The number of dimensions DimSizes must contain ndim elements
     * @param Driver The function object that determines whether we do forward or inverse
     */
    template<class WaveletDriverType>
    void WaveletTransformImp(double *InArray, const size_t *DimSizes,
        size_t ndim, WaveletDriverType Driver)
      {
        const size_t nElements = std::accumulate(DimSizes, DimSizes + ndim, 1,
            std::multiplies<size_t>());
        jif3D::rvec WorkVector(nElements);
        size_t Stride = 1;
        size_t preStride = 1;
        for (size_t i = 0; i < ndim; ++i)
          {
            size_t DimSize = DimSizes[i];
            Stride = DimSize * preStride;
            if (DimSize > 4)
              {
                for (size_t j = 0; j < nElements; j += Stride)
                  {
                    for (size_t k = 0; k < preStride; ++k)
                      {
                        for (size_t l = j + k, m = 0; m < DimSize; ++m, l
                            += preStride)
                          {
                            WorkVector(m) = InArray[l];
                          }

                        Driver(WorkVector, DimSize);

                        for (size_t l = j + k, m = 0; m < DimSize; ++m, l
                            += preStride)
                          {
                            InArray[l] = WorkVector(m);
                          }
                      }
                  }
              }
            preStride = Stride;
          }
      }

    /*! Perform a multidimensional wavelet transform on the input data. InArray is a one dimensional array
     * of DimSizes[0]*DimSizes[1]* ... *DimSizes[ndim-1] elements. We assume that the last dimension varies
     * fastest (c storage order).
     * @param InArray The array with the input values, contains the output afterwards
     * @param DimSizes The size of each dimension, InArray must contain the product of DimSizes elements
     * @param ndim The number of dimensions DimSizes must contain ndim elements
     */
    void WaveletTransform(double *InArray, const size_t *DimSizes, size_t ndim)
      {
        //we make sure that each dimension has a size that is a power of two
        for (size_t i = 0; i < ndim; ++i)
          {
            if (!IsPowerOfTwo(DimSizes[i]))
              throw jif3D::FatalException(
                  "Array dimension is not a power of two.");
          }
        WaveletTransformImp(InArray, DimSizes, ndim, ForwardWaveletDriver());
      }

    void InvWaveletTransform(double *InArray, const size_t *DimSizes,
        size_t ndim)
      {
        //we make sure that each dimension has a size that is a power of two
        for (size_t i = 0; i < ndim; ++i)
          {
            if (!IsPowerOfTwo(DimSizes[i]))
              throw jif3D::FatalException(
                  "Array dimension is not a power of two.");
          }
        WaveletTransformImp(InArray, DimSizes, ndim, InverseWaveletDriver());
      }

  }
