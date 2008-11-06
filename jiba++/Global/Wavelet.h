//============================================================================
// Name        : Wavelet.h
// Author      : Sep 26, 2008
// Version     :
// Copyright   : 2008, mmoorkamp
//============================================================================

#include "VecMat.h"
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
  /* @} */
  }
