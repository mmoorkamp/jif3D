/*
 * X3DFreqFunctions.h
 *
 *  Created on: Mar 17, 2014
 *      Author: mmoorkamp
 */

#ifndef X3DFREQFUNCTIONS_H_
#define X3DFREQFUNCTIONS_H_

#ifdef HAVEHPX
#include <hpx/hpx.hpp>
#endif

#include <string>
#include <vector>

#include "../Global/Jif3DGlobal.h"
#include "../Global/Serialization.h"
#include "../Global/VecMat.h"
#include "../Global/VectorTransform.h"
#include "ReadWriteX3D.h"
#include "X3DModel.h"
#include "X3DFieldCalculator.h"
#include "MTData.h"
#include "TipperData.h"
#include "X3DTypes.h"


ForwardResult CalculateFrequency(const ForwardInfo &Info, const jif3D::MTData &Data,
    const jif3D::X3DFieldCalculator &Calc);

GradResult LQDerivativeFreq(const ForwardInfo &Info, const  jif3D::MTData &Data, const GradInfo &GI,
    const jif3D::X3DFieldCalculator &Calc);

GradResult TipperDerivativeFreq(const ForwardInfo &Info, const jif3D::TipperData &Data, const jif3D::rvec &Misfit,
    const jif3D::X3DFieldCalculator &Calc);

#ifdef HAVEHPX
HPX_DEFINE_PLAIN_ACTION(CalculateFrequency, CalculateFrequency_action);
HPX_DEFINE_PLAIN_ACTION(LQDerivativeFreq, LQDerivativeFreq_action);
#endif

#endif /* X3DFREQFUNCTIONS_H_ */
