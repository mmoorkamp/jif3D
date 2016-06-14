/*
 * CalcRealization.h
 *
 *  Created on: 19 Aug 2014
 *      Author: mm489
 */

#ifndef CALCREALIZATION_H_
#define CALCREALIZATION_H_

#ifdef HAVEHPX
#include <hpx/hpx.hpp>
#include <hpx/include/actions.hpp>
#endif
#include "../Global/Serialization.h"
#include "../Global/VecMat.h"
#include "../Global/Jif3DGlobal.h"
#include "realinfo.h"

J3DEXPORT std::vector<double> CalcRealization(realinfo Info);

#ifdef HAVEHPX
HPX_DEFINE_PLAIN_ACTION(CalcRealization, Calc_action);
#endif

#endif /* CALCREALIZATION_H_ */
