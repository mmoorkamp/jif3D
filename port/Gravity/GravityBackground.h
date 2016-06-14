//============================================================================
// Name        : GravityBackground.h
// Author      : Jun 18, 2009
// Version     : 
// Copyright   : 2009, mmoorkamp
//============================================================================

#ifndef GRAVITYBACKGROUND_H_
#define GRAVITYBACKGROUND_H_

namespace jif3D
  {
    //! Calculate the effect of a 1D layered background on scalar gravity data
    rvec CalcScalarBackground(const size_t measindex, const double xwidth,
        const double ywidth, const double zwidth, const ThreeDGravityModel &Model,
        rmat &Sensitivities);
    //! Calculate the effect of a 1D layered background on full tensor gravity (FTG) data
    rvec CalcTensorBackground(const size_t measindex, const double xwidth,
        const double ywidth, const double zwidth, const ThreeDGravityModel &Model,
        rmat &Sensitivities);
  }

#endif /* GRAVITYBACKGROUND_H_ */
