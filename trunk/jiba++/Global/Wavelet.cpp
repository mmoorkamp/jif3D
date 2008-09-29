//============================================================================
// Name        : Wavelet.cpp
// Author      : Sep 26, 2008
// Version     :
// Copyright   : 2008, mmoorkamp
//============================================================================

#include "Wavelet.h"
namespace jiba
  {

    static const double c0 = (1.0 + sqrt(3.0))/(4.0*sqrt(2.0));
    static const double c1 = (3.0+sqrt(3.0))/(4.0*sqrt(2));
    static const double c2 = (3.0-sqrt(3.0))/(4.0*sqrt(2.0));
    static const double c3 = (1.0-sqrt(3.0))/(4.0*sqrt(2.0));

    void Daub4(jiba::rvec &Invec, const size_t maxindex)
      {
        jiba::rvec Temp(Invec);
        const size_t length = maxindex;
        const size_t offset = length / 2;
        for (size_t i = 0, j = 0; j < length - 3; j += 2, ++i)
          {
            Temp(i) = c0 * Invec(j) + c1 * Invec(j+1) + c2 * Invec(j+2) + c3 * Invec(j+3);
            Temp(i + offset) = c3 * Invec(j) - c2 * Invec(j + 1) + c1 * Invec(j
                + 2) - c0 * Invec(j + 3);
          }
        Temp(length / 2 - 1) = c0 * Invec(length - 2) + c1 * Invec(length - 1)
            + c2 * Invec(0) + c3 * Invec(1);
        Temp(length - 1) = c3 * Invec(length - 2) - c2 * Invec(length - 1) + c1
            * Invec(0) - c0 * Invec(1);
        std::copy(Temp.begin(), Temp.begin() + maxindex, Invec.begin());
      }

    void InvDaub4(jiba::rvec &Invec, const size_t maxindex)
      {
        jiba::rvec Temp(Invec);
        const size_t length = maxindex;
        const size_t offset = length / 2;
        Temp(0) = c2 * Invec(length / 2 - 1) + c1 * Invec(length - 1) + c0
            * Invec(0) + c3 * Invec(length / 2);
        Temp(1) = c3 * Invec(length / 2 - 1) - c0 * Invec(length - 1) + c1
            * Invec(0) - c2 * Invec(length / 2);
        for (size_t i = 0, j = 2; i < offset - 1; ++i, j += 2)
          {
            Temp(j) = c2 * Invec(i) + c1 * Invec(i+offset) + c0
            * Invec(i+1) + c3 * Invec(i+offset+1);
            Temp(j + 1) = c3 * Invec(i) - c0 * Invec(i + offset) + c1 * Invec(i
                + 1) - c2 * Invec(i + offset + 1);
          }
        std::copy(Temp.begin(), Temp.begin() + maxindex, Invec.begin());
      }
    void WaveletTransform(jiba::rvec &Invec)
      {
        const size_t length = Invec.size();
        for (size_t i = length; i >= 4; i /= 2)
          {
            Daub4(Invec, i);
          }
      }
    void InvWaveletTransform(jiba::rvec &Invec)
      {
        const size_t length = Invec.size();
        for (size_t i = 4; i <= length; i *= 2)
          {
            InvDaub4(Invec, i);
          }
      }
  }
