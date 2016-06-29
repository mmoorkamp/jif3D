//============================================================================
// Name        : Tarits2DMT.h
// Author      : Max Moorkamp
// Version     :
// Copyright   : 2008, MM
//============================================================================

#ifndef TARITS2DMT_H_
#define TARITS2DMT_H_

#include "../Global/Jif3DGlobal.h"

extern "C"
  {
    /** \addtogroup mtmodelling Forward modelling of magnetotelluric data */
    /* @{ */

    //! Calculate the MT H-Polarization response for a 2D resistivity model, this is an interface to Tarits's Fortran routines
    /*!
     * @param T The period in s
     * @param NX The number of cells in x-direction
     * @param NZ The number of cells in z-direction
     * @param DXJ Vector of cell sizes in m in x-direction
     * @param DZI Vector of cell sizes in m in z-direction
     * @param RO Vector of size nx*nz that contains the resistivity of each cell in Ohmm, c-style storage
     * @param HY_REAL The real part of the magnetic field in each cell (output)
     * @param HY_IMAG The imaginary part of the magnetic field in each cell (output)
     * @param EX_REAL The real part of the horizontal electric field in each cell (output)
     * @param EX_IMAG The imaginary part of the horizontal electric field in each cell (output)
     * @param EZ_REAL The real part of the vertical electric field in each cell (output)
     * @param EZ_IMAG The imaginary part of the vertical electric field in each cell (output)
     */
    J3DEXPORT void hpol_(const double* T, const long* NX, const long* NZ,
        const double* DXJ, const double* DZI, const double* RO,
        double* HY_REAL, double* HY_IMAG, double* EX_REAL, double* EX_IMAG,
        double* EZ_REAL, double* EZ_IMAG);

    //! Calculate the MT E-Polarization response for a 2D resistivity model, this is an interface to Tarits's Fortran routines
    /*!
     *
     * @param T The period in s
     * @param NX The number of cells in x-direction
     * @param NZ The number of cells in z-direction
     * @param IONOS The number of ionospheric layers
     * @param IATM The number of atmospheric layers
     * @param RIONOS The resistivity of the ionosphere
     * @param DXJ Vector of cell sizes in m in x-direction
     * @param DZI Vector of cell sizes in m in z-direction. This includes the air and ionospheric layers !
     * @param RO Vector of size nx*nz that contains the resistivity of each cell in Ohmm, c-style storage
     * @param HX_REAL The real part of the magnetic field in each cell (output)
     * @param HX_IMAG The imaginary part of the magnetic field in each cell (output)
     * @param EY_REAL The real part of the horizontal electric field in each cell (output)
     * @param EY_IMAG The imaginary part of the horizontal electric field in each cell (output)
     */
    J3DEXPORT void epol_(const double* T, const long* NX, const long* NZ,
        const long* IONOS, const long* IATM, const double* DXJ,
        const double* DZI, const double* RO, const double* RIONOS,
        double* HX_REAL, double* HX_IMAG, double* EY_REAL, double* EY_IMAG);
  /* @} */
  }

#endif /*TARITS2DMT_H_*/

