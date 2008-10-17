//============================================================================
// Name        : Wavelet.h
// Author      : Sep 26, 2008
// Version     :
// Copyright   : 2008, mmoorkamp
//============================================================================

#include "VecMat.h"
namespace jiba
  {
    /** \addtogroup util General utility routines */
    /* @{ */
    void Daub4(jiba::rvec &Invec, const size_t maxindex);
    void InvDaub4(jiba::rvec &Invec, const size_t maxindex);
    void WaveletTransform(jiba::rvec &Invec);
    void InvWaveletTransform(jiba::rvec &Invec);
  /* @} */
  }
