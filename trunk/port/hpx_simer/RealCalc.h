/*
 * RealCalc.h
 *
 *  Created on: 19 Aug 2014
 *      Author: mm489
 */

#ifndef REALCALC_H_
#define REALCALC_H_

#include "realinfo.h"
#include "../Global/Jif3DGlobal.h"

class J3DEXPORT RealCalc
  {
  realinfo Info;
public:
  void SetInfo(realinfo I)
    {
      Info = I;
    }
  RealCalc();
  virtual ~RealCalc();
  };

#endif /* REALCALC_H_ */
