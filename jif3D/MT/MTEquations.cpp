//============================================================================
// Name        : MTEquations.cpp
// Author      : Jul 7, 2009
// Version     : 
// Copyright   : 2009, mmoorkamp
//============================================================================

#include "MTEquations.h"
#include "../Global/convert.h"
#include "../Global/FatalException.h"

namespace jif3D
  {
    /*! Transform the field spectra to impedance values in S.I. units.
     * @param Ex1 The electric field in x-direction for the first source polarization
     * @param Ex2 The electric field in x-direction for the second source polarization
     * @param Ey1 The electric field in y-direction for the first source polarization
     * @param Ey2 The electric field in y-direction for the second source polarization
     * @param Hx1 The magnetic field in x-direction for the first source polarization
     * @param Hx2 The magnetic field in x-direction for the second source polarization
     * @param Hy1 The magnetic field in y-direction for the first source polarization
     * @param Hy2 The magnetic field in y-direction for the second source polarization
     * @param Zxx The xx element of the impedance in S.I. units
     * @param Zxy The xy element of the impedance in S.I. units
     * @param Zyx The yx element of the impedance in S.I. units
     * @param Zyy The yy element of the impedance in S.I. units
     */
    void FieldsToImpedance(const std::complex<double> &Ex1,
        const std::complex<double> &Ex2, const std::complex<double> &Ey1,
        const std::complex<double> &Ey2, const std::complex<double> &Hx1,
        const std::complex<double> &Hx2, const std::complex<double> &Hy1,
        const std::complex<double> &Hy2, std::complex<double> &Zxx,
        std::complex<double> &Zxy, std::complex<double> &Zyx, std::complex<double> &Zyy)
      {
        const std::complex<double> magdet(Hx1 * Hy2 - Hy1 * Hx2);
        Zxx = (Ex1 * Hy2 - Hy1 * Ex2) / magdet;
        Zxy = (Ex2 * Hx1 - Hx2 * Ex1) / magdet;
        Zyx = (Ey1 * Hy2 - Hy1 * Ey2) / magdet;
        Zyy = (Ey2 * Hx1 - Hx2 * Ey1) / magdet;
      }

    /*! For a given frequency and conductivity, calculate the impedance
     * for a half-space with that conductivity.
     * @param frequency The measurement frequency in Hz
     * @param conductivity The conductivity of the half-space in S/m
     * @return The impedance of the half-space in S.I. units, i.e. Ohm
     */
    std::complex<double> ImpedanceHalfspace(const double frequency,
        const double conductivity)
      {
        std::complex<double> omegamu = std::complex<double>(0.0, twopimu * frequency);
        return sqrt(omegamu / conductivity);
      }

    /*! Rotate the impedance elements by the given angle in radian, the input variables
     * will contain the rotated values;
     * @param angle Rotation angle clockwise from north in radian
     * @param Zxx Zxx element of the impedance tensor
     * @param Zxy Zxy element of the impedance tensor
     * @param Zyx Zyx element of the impedance tensor
     * @param Zyy Zyy element of the impedance tensor
     */
    void RotateImpedance(const double angle, std::complex<double> & Zxx,
        std::complex<double> &Zxy, std::complex<double> &Zyx, std::complex<double> &Zyy)
      {
        //we need the old impedance elements in all 4 equations
        //so we create some temporary variable
        std::complex<double> newxx, newxy, newyx, newyy;

        //precompute trigonometric expressions used several times below
        const double ca2 = pow(cos(angle), 2);
        const double sa2 = pow(sin(angle), 2);
        const double casa = sin(angle) * cos(angle);
        //do the rotation R Z R^T
        newxx = Zxx * ca2 - (Zxy + Zyx) * casa + Zyy * sa2;
        newxy = Zxy * ca2 + (Zxx - Zyy) * casa - Zyx * sa2;
        newyx = Zyx * ca2 + (Zxx - Zyy) * casa - Zxy * sa2;
        newyy = Zyy * ca2 + (Zxy + Zyx) * casa + Zxx * sa2;
        //assign temporary values to impedance elements
        Zxx = newxx;
        Zxy = newxy;
        Zyx = newyx;
        Zyy = newyy;
      }

    jif3D::rvec RotateImpedanceVector(const double angle, const jif3D::rvec &Impedance)
      {
        jif3D::rvec Result(Impedance.size());
        const size_t nelem = 8;
        if (Impedance.size() % nelem != 0)
          {
            throw jif3D::FatalException(
                "Size of impedance vector: " + stringify(Impedance.size())
                    + " is not a multiple of 8!", __FILE__, __LINE__);
          }
        const size_t ntensors = Impedance.size() / nelem;
        for (size_t i = 0; i < ntensors; ++i)
          {
            std::complex<double> Zxx(Impedance(i * nelem), Impedance(i * nelem + 1));
            std::complex<double> Zxy(Impedance(i * nelem + 2), Impedance(i * nelem + 3));
            std::complex<double> Zyx(Impedance(i * nelem + 4), Impedance(i * nelem + 5));
            std::complex<double> Zyy(Impedance(i * nelem + 6), Impedance(i * nelem + 7));
            RotateImpedance(angle, Zxx, Zxy, Zyx, Zyy);
            Result(i * nelem) = Zxx.real();
            Result(i * nelem + 1) = Zxx.imag();
            Result(i * nelem + 2) = Zxy.real();
            Result(i * nelem + 3) = Zxy.imag();
            Result(i * nelem + 4) = Zyx.real();
            Result(i * nelem + 5) = Zyx.imag();
            Result(i * nelem + 6) = Zyy.real();
            Result(i * nelem + 7) = Zyy.imag();
          }
        return Result;
      }

  }
