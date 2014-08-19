/*
 * CalcRealization.h
 *
 *  Created on: 19 Aug 2014
 *      Author: mm489
 */

#ifndef CALCREALIZATION_H_
#define CALCREALIZATION_H_

#ifdef HAVEHPX
#include <hpx/config.hpp>
#include <hpx/include/actions.hpp>
#endif
#include "../Global/VecMat.h"
#include "realinfo.h"


jif3D::rvec CalcRealization(realinfo Info);

#ifdef HAVEHPX
HPX_DEFINE_PLAIN_ACTION(CalcRealization, Calc_action);
#endif




#endif /* CALCREALIZATION_H_ */
