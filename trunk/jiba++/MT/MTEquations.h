//============================================================================
// Name        : MTEquations.h
// Author      : Jul 7, 2009
// Version     : 
// Copyright   : 2009, mmoorkamp
//============================================================================

#include <complex>
namespace jiba
  {
    const double mag_mu = 4e-7 * M_PI;

    std::complex<double> ImpedanceHalfspace(const double frequency,
        const double conductivity);

    void FieldsToImpedance(const std::complex<double> &Ex1, const std::complex<
        double> &Ex2, const std::complex<double> &Ey1, const std::complex<
        double> &Ey2, const std::complex<double> &Hx1, const std::complex<
        double> &Hx2, const std::complex<double> &Hy1, const std::complex<
        double> &Hy2, std::complex<double> &Zxx, std::complex<double> &Zxy,
        std::complex<double> &Zyx, std::complex<double> &Zyy);

  }

