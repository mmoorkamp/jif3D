//============================================================================
// Name        : Wavelet.h
// Author      : Sep 26, 2008
// Version     :
// Copyright   : 2008, mmoorkamp
//============================================================================

#include "VecMat.h"
#include <boost/multi_array.hpp>

/*! \file Wavelet.h
 * Simple implementations of the discrete wavelet transform and its inverse based on Daubechies-4 wavelets.
 * The basic algorithm is described in Numerical Recipes in C, pp. 595-598
 */

namespace jiba
  {
    /** \addtogroup util General utility routines */
    /* @{ */
    //! Apply the Daubechies 4 coefficient wavelet filter to a vector
    void Daub4(jiba::rvec &Invec, const size_t maxindex);
    //! Apply the inverse Daubechies 4 coefficient wavelet filter to a vector
    void InvDaub4(jiba::rvec &Invec, const size_t maxindex);
    //! Perform a forward wavelet transform on the input vector
    void WaveletTransform(jiba::rvec &Invec);
    //! Perform an inverse wavelet transform on the input vector
    void InvWaveletTransform(jiba::rvec &Invec);
    //! Multidimensional wavelet transform
    template<typename MultiArrayType>
    void WaveletTransform(MultiArrayType &InArray);
    template<typename MultiArrayType>
    void InvWaveletTransform(MultiArrayType &InArray);
    /* @} */

    template<typename MultiArrayType>
    void WaveletTransform(MultiArrayType &InArray)
      {
        const size_t nDimensions = InArray.num_dimensions();
        const size_t nElements = std::accumulate(InArray.shape(),
            InArray.shape() + nDimensions, 1,std::multiplies<size_t>());
        jiba::rvec WorkVector(nElements);
        size_t Stride = 1;
        size_t preStride = 1;
        for (size_t i = 0; i < nDimensions; ++i)
          {
            size_t DimSize = InArray.shape()[i];
            Stride = DimSize * preStride;
            if (InArray.shape()[i] > 4)
              {
                for (size_t j = 0; j < nElements; j += Stride)
                  {
                    for (size_t k = 0; k < preStride; ++k)
                      {
                        for (size_t l = j + k, m = 0; m < DimSize; ++m, l
                            += preStride)
                          {
                            WorkVector( m) = *(InArray.origin() + l);
                          }
                        for (size_t n = DimSize; n >= 4; n >>= 1)
                          {
                            Daub4(WorkVector, n);
                          }
                        for (size_t l = j + k, m = 0; m < DimSize; ++m, l
                            += preStride)
                          {
                            *(InArray.origin() + l) = WorkVector(m);
                          }
                      }
                  }
              }
            preStride = Stride;
          }
      }

    template<typename MultiArrayType>
        void InvWaveletTransform(MultiArrayType &InArray)
          {
            const size_t nDimensions = InArray.num_dimensions();
            const size_t nElements = std::accumulate(InArray.shape(),
                InArray.shape() + nDimensions, 1,std::multiplies<size_t>());
            jiba::rvec WorkVector(nElements);
            size_t Stride = 1;
            size_t preStride = 1;
            for (size_t i = 0; i < nDimensions; ++i)
              {
                size_t DimSize = InArray.shape()[i];
                Stride = DimSize * preStride;
                if (InArray.shape()[i] > 4)
                  {
                    for (size_t j = 0; j < nElements; j += Stride)
                      {
                        for (size_t k = 0; k < preStride; ++k)
                          {
                            for (size_t l = j + k, m = 0; m < DimSize; ++m, l
                                += preStride)
                              {
                                WorkVector( m) = *(InArray.origin() + l);
                              }
                            for (size_t n = 4; n <= DimSize; n <<= 1)
                              {
                                InvDaub4(WorkVector, n);
                              }
                            for (size_t l = j + k, m = 0; m < DimSize; ++m, l
                                += preStride)
                              {
                                *(InArray.origin() + l) = WorkVector(m);
                              }
                          }
                      }
                  }
                preStride = Stride;
              }
          }
  }
