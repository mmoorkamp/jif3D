#ifndef TARITS2DMT_H_
#define TARITS2DMT_H_

/*Fortran Routines (modified from Tartis)*/
/*in files B_hpol.for and B_epol.for*/

//! Calculate the MT H-Polarization response for a 2D resistivity model
extern "C" void HPOL(double* T, long* NX, long* NZ, double* DXJ, double* DZI,
    double* RO, double* HY_REAL, double* HY_IMAG, double* EX_REAL,
    double* EX_IMAG, double* EZ_REAL, double* EZ_IMAG);

//! Calculate the MT E-Polarization response for a 2D resistivity model
extern "C" void epol_(double* T, int* NX, int* NZ, int* IONOS, int* IATM,
    double* DXJ, double* DZI, double* RO, double* RIONOS, double* HX_REAL,
    double* HX_IMAG, double* EY_REAL, double* EY_IMAG);

#endif /*TARITS2DMT_H_*/
