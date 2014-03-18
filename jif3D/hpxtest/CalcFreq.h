//============================================================================
// Name        : CalcFreq.h
// Author      : 23 Feb 2013
// Version     : 
// Copyright   : 2013, mm489
//============================================================================


#ifndef CALCFREQ_H_
#define CALCFREQ_H_

#include <hpx/config.hpp>
#include <hpx/include/actions.hpp>
#include "../MT/X3DModel.h"
#include "../Global/VecMat.h"


jif3D::rvec HPXCalculateFrequency(const jif3D::X3DModel &Model, size_t freqindex,
    std::string TempDirName, std::string X3DName);

HPX_DEFINE_PLAIN_ACTION(HPXCalculateFrequency, CalculateFrequency_action);
//HPX_ACTION_USES_SNAPPY_COMPRESSION(CalculateFrequency_action);

#endif /* CALCFREQ_H_ */
