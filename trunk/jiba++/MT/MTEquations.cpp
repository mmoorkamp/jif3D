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
        /*if (abs(Hx1) == 0 && abs(Hy2) == 0)
         {
         Zxx = std::complex<double>();
         Zyy = std::complex<double>();
         Zxy = Ex1 / Hy1;
         Zyx = Ey2 / Hx2;
         }*/
        const std::complex<double> magdet(Hx1 * Hy2 - Hy1 * Hx2);
        Zxx = (Ex1 * Hy2 - Hy1 * Ex2) / magdet;
        Zxy = (Ex2 * Hx1 - Hx2 * Ex1) / magdet;
        Zyx = (Ey1 * Hy2 - Hy1 * Ey2) / magdet;
        Zyy = (Ey2 * Hx1 - Hx2 * Ey1) / magdet;
      }


    //! Calculate the impedance of a half-space in S.I. units (Ohm)
    /*! For a given frequency and conductivity, calculate the impedance
     * for a half-space with that conductivity.
     * @param frequency The measurement frequency in Hz
     * @param conductivity The conductivity of the half-space in S/m
     * @return The impedance of the half-space in S.I. units, i.e. Ohm
     */
    std::complex<double> ImpedanceHalfspace(const double frequency,
        const double conductivity)
      {
        std::complex<double> omegamu = std::complex<double>(0.0, 2.0 * M_PI
            * mag_mu * frequency);
        return sqrt(omegamu / conductivity);
      }
  }
