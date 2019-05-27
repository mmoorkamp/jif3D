//============================================================================
// Name        : GravityBackground.h
// Author      : Jun 18, 2009
// Version     :
// Copyright   : 2009, mmoorkamp
//============================================================================

#ifndef GRAVITYBACKGROUND_H_
#define GRAVITYBACKGROUND_H_

#include "../Global/Jif3DGlobal.h"
#include "ThreeDGravityModel.h"
#include "ScalarGravityData.h"
#include "TensorGravityData.h"

namespace jif3D
  {
    //! Calculate the effect of a 1D layered background on scalar gravity data
    J3DEXPORT rvec CalcScalarBackground(const size_t measindex, const double xwidth,
        const double ywidth, const double zwidth, const ThreeDGravityModel &Model, const ScalarGravityData &Data,
        rmat &Sensitivities);
    //! Calculate the effect of a 1D layered background on full tensor gravity (FTG) data
    J3DEXPORT rvec CalcTensorBackground(const size_t measindex, const double xwidth,
        const double ywidth, const double zwidth, const ThreeDGravityModel &Model,
        const TensorGravityData &Data,
        rmat &Sensitivities);
  }

#endif /* GRAVITYBACKGROUND_H_ */
