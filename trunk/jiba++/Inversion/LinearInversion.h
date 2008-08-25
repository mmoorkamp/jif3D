//============================================================================
// Name        : gravgrid.cpp
// Author      : Max Moorkamp
// Version     :
// Copyright   : 2008, MM
//============================================================================


#ifndef LINEARINVERSION_H_
#define LINEARINVERSION_H_

#include "../Global/VecMat.h"

namespace jiba
  {
    void DataSpaceInversion(const rmat &Sensitivities, const rvec &Data,
        const rvec &WeightVector, const double evalthresh, rvec &InvModel);
  }
#endif /* LINEARINVERSION_H_ */
