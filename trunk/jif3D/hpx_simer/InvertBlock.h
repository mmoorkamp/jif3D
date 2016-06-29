//============================================================================
// Name        : InvertBlock.h
// Author      : 19 Dec 2014
// Version     :
// Copyright   : 2014, mm489
//============================================================================

#ifndef HPX_SIMER_INVERTBLOCK_H_
#define HPX_SIMER_INVERTBLOCK_H_


#include <vector>
#include <string>
#include "../Global/Jif3DGlobal.h"
#include "../Global/VecMat.h"
#include "../MT/X3DModel.h"
#include "realinfo.h"

J3DEXPORT double InvertBlock(jif3D::X3DModel Model, jif3D::rvec Data,  std::string tempdir, std::string x3dname);

#endif /* HPX_SIMER_INVERTBLOCK_H_ */
