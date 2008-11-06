//============================================================================
// Name        : Wavelet.cpp
// Author      : Sep 26, 2008
// Version     :
// Copyright   : 2008, mmoorkamp
//============================================================================

#include "Wavelet.h"
namespace jiba
  {
    //calcualte the coefficient for the wavelet filter
    static const double c0 = (1.0 + sqrt(3.0)) / (4.0 * sqrt(2.0));
    static const double c1 = (3.0 + sqrt(3.0)) / (4.0 * sqrt(2));
    static const double c2 = (3.0 - sqrt(3.0)) / (4.0 * sqrt(2.0));
    static const double c3 = (1.0 - sqrt(3.0)) / (4.0 * sqrt(2.0));
    /*! Apply a 4 coefficient Daubechies wavelet filter to the input vector.
     * @param Invec The vector with the input data, will contain the result of the transform. Size must be a power of 2
     * @param maxindex The maximum index to which to apply the filter, must be a power of 2
     */
    void Daub4(jiba::rvec &Invec, const size_t maxindex)
      {
        //we perform the operations on a temporary copy
        jiba::rvec Temp(Invec);
        const size_t length = maxindex;
        const size_t offset = length / 2;
        //the current value depends on the next three values, so we have to adjust the end of the loop
        //also note the different stride for j and i
        for (size_t i = 0, j = 0; j < length - 3; j += 2, ++i)
          {
            Temp( i) = c0 * Invec(j) + c1 * Invec(j + 1) + c2 * Invec(j + 2)
                + c3 * Invec(j + 3);
            Temp(i + offset) = c3 * Invec(j) - c2 * Invec(j + 1) + c1 * Invec(j
                + 2) - c0 * Invec(j + 3);
          }
        //we have to handle the last values as a special case
        Temp(length / 2 - 1) = c0 * Invec(length - 2) + c1 * Invec(length - 1)
            + c2 * Invec(0) + c3 * Invec(1);
        Temp(length - 1) = c3 * Invec(length - 2) - c2 * Invec(length - 1) + c1
            * Invec(0) - c0 * Invec(1);
        //copy the result back to Invec
        std::copy(Temp.begin(), Temp.begin() + maxindex, Invec.begin());
      }
    /*! Apply a 4 coefficient inverse Daubechies wavelet filter to the input vector.
     * @param Invec The vector with the input data, will contain the result of the transform. Size must be a power of 2
     * @param maxindex The maximum index to which to apply the filter, must be a power of 2
     */
    void InvDaub4(jiba::rvec &Invec, const size_t maxindex)
      {
        //this is the same strategy as the forward
        //only the equations and the order are slightly different
        jiba::rvec Temp(Invec);
        const size_t length = maxindex;
        const size_t offset = length / 2;
        Temp(0) = c2 * Invec(length / 2 - 1) + c1 * Invec(length - 1) + c0
            * Invec(0) + c3 * Invec(length / 2);
        Temp(1) = c3 * Invec(length / 2 - 1) - c0 * Invec(length - 1) + c1
            * Invec(0) - c2 * Invec(length / 2);
        for (size_t i = 0, j = 2; i < offset - 1; ++i, j += 2)
          {
            Temp( j) = c2 * Invec(i) + c1 * Invec(i + offset) + c0 * Invec(i
                + 1) + c3 * Invec(i + offset + 1);
            Temp(j + 1) = c3 * Invec(i) - c0 * Invec(i + offset) + c1 * Invec(i
                + 1) - c2 * Invec(i + offset + 1);
          }
        std::copy(Temp.begin(), Temp.begin() + maxindex, Invec.begin());
      }
    /*! Perform a forward wavelet transform on the input vector
     * @param Invec The input vector, will contain the result, size must be power of 2
     */
    void WaveletTransform(jiba::rvec &Invec)
      {
        const size_t length = Invec.size();
        for (size_t i = length; i >= 4; i /= 2)
          {
            Daub4(Invec, i);
          }
      }
    /*! Perform an inverse wavelet transform on the input vector
     * @param Invec The input vector, will contain the result, size must be power of 2
     */
    void InvWaveletTransform(jiba::rvec &Invec)
      {
        const size_t length = Invec.size();
        for (size_t i = 4; i <= length; i *= 2)
          {
            InvDaub4(Invec, i);
          }
      }
  }
