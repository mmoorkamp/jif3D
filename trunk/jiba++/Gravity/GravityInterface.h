//============================================================================
// Name        : GravityInterface.h
// Author      : Max Moorkamp
// Version     :
// Copyright   : 2008, MM
//============================================================================

#ifndef GRAVITYINTERFACE_H_
#define GRAVITYINTERFACE_H_

/*! \file Here we declare the interface to the 3D Gravity code for use with R
 */

extern "C"
void  CalcScalarForward(const double *XSizes, const unsigned int nx, const double *YSizes, const unsigned int ny,
      const double *ZSizes, const unsigned int nz, const double *Densities,
      const double *XMeasPos, const double *YMeasPos, const double *ZMeasPos, const unsigned int nmeas,
      double *GravAcceleration);

extern "C"
void  CalcTensorForward();
#endif /* GRAVITYINTERFACE_H_ */
