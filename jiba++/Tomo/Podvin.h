//============================================================================
// Name        : Podvin.h
// Author      : Sep 24, 2008
// Version     :
// Copyright   : 2008, mmoorkamp
//============================================================================


#ifndef PODVIN_H_
#define PODVIN_H_

extern "C"
  {
    int time_3d(float *HS, float *T, int NX, int NY, int NZ, float XS,
        float YS, float ZS, float HS_EPS_INIT, int MSG);
  }
#endif /* PODVIN_H_ */
