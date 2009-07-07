//============================================================================
// Name        : MTEquations.cpp
// Author      : Jul 7, 2009
// Version     : 
// Copyright   : 2009, mmoorkamp
//============================================================================


#include "MTEquations.h"

namespace jiba
  {

    void FieldsToImpedance(const std::complex<double> &Ex1, const std::complex<
        double> &Ex2, const std::complex<double> &Ey1, const std::complex<
        double> &Ey2, const std::complex<double> &Hx1, const std::complex<
        double> &Hx2, const std::complex<double> &Hy1, const std::complex<
        double> &Hy2, std::complex<double> &Zxx, std::complex<double> &Zxy,
        std::complex<double> &Zyx, std::complex<double> &Zyy)
      {
        //1D case
        if (abs(Hx1) == 0 && abs(Hy2) == 0 )
          {
            Zxx = std::complex<double>();
            Zyy = std::complex<double>();
            Zxy = Ex1/Hy1;
            Zyx = Ey2/Hx2;
          }
        const std::complex<double> magdet(Hx1 * Hy2 - Hy1 * Hx2);
        Zxx = (Ex1 * Hy2 - Hy1 * Ex2) / magdet;
        Zxy = (Ex2 * Hx1 - Hx2 * Ex1) / magdet;
        Zyx = (Ey1 * Hy2 - Hy1 * Ey2) / magdet;
        Zyy = (Ey1 * Hx1 - Hx2 * Ey1) / magdet;
      }

    std::complex<double> ImpedanceHalfspace(const double frequency,
        const double conductivity)
      {
        std::complex<double> omegamu = std::complex<double>(0.0, 1.0) * 8e-7
            * M_PI * M_PI * frequency;
        return sqrt(omegamu / conductivity);
      }
  }
