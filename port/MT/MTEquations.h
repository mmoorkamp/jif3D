//============================================================================
// Name        : MTEquations.h
// Author      : Jul 7, 2009
// Version     :
// Copyright   : 2009, mmoorkamp
//============================================================================

#include <complex>

#include <boost/math/constants/constants.hpp>

#include "../Global/Jif3DGlobal.h"
#include "../Global/VecMat.h"

namespace jif3D
  {
    /** \addtogroup mtmodelling Forward modelling of magnetotelluric data */
    /* @{ */
    //! The value of magnetic permeability \f$ \mu = 4 \cdot 10^{-7} \pi\f$
    const double mag_mu = 4e-7 * boost::math::constants::pi<double>();
    //! Two pi times the magnetic permeability, needed for apparent resistivity calculation below
    const double twopimu = 2.0 * boost::math::constants::pi<double>() * mag_mu;
    //! Calculate the impedance of a half-space in S.I. units (Ohm)
    J3DEXPORT std::complex<double> ImpedanceHalfspace(const double frequency,
        const double conductivity);
    //! From the spectra of the electric and magnetic fields, calculate the impedance tensor elements
    J3DEXPORT void FieldsToImpedance(const std::complex<double> &Ex1,
        const std::complex<double> &Ex2, const std::complex<double> &Ey1,
        const std::complex<double> &Ey2, const std::complex<double> &Hx1,
        const std::complex<double> &Hx2, const std::complex<double> &Hy1,
        const std::complex<double> &Hy2, std::complex<double> &Zxx,
        std::complex<double> &Zxy, std::complex<double> &Zyx, std::complex<double> &Zyy);
    //! Rotate the impedance tensor elements by the given angle in radian
    J3DEXPORT void RotateImpedance(const double angle, std::complex<double> & Zxx,
        std::complex<double> &Zxy, std::complex<double> &Zyx, std::complex<double> &Zyy);
    //! Rotate several impedance values stored in a vector as we use for inversion
    J3DEXPORT jif3D::rvec RotateImpedanceVector(const double angle, const jif3D::rvec &Impedance);
    //! Given a complex impedance value in S.I. units, and the Frequency in Hz calculate the corresponding apparent resistivity
    inline double AppRes(const std::complex<double> &Z, const double Frequency)
      {
        return 1.0 / (twopimu * Frequency) * std::norm(Z);
      }
    //! Return the phase of the magnetotelluric impedance in degree
    inline double ImpedancePhase(const std::complex<double> &Z)
      {
        return 180.0 / boost::math::constants::pi<double>() * std::atan2(Z.imag(), Z.real());
      }
  /* @} */
  }

