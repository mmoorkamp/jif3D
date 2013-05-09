//============================================================================
// Name        : CalcFreq.h
// Author      : 23 Feb 2013
// Version     : 
// Copyright   : 2013, mm489
//============================================================================


#ifndef CALCFREQ_H_
#define CALCFREQ_H_

#include <hpx/include/lcos.hpp>
#include <hpx/include/iostreams.hpp>
#include <hpx/include/actions.hpp>
#include <hpx/include/util.hpp>
#include "../MT/X3DModel.h"
#include "../Global/VecMat.h"


jiba::rvec HPXCalculateFrequency(const jiba::X3DModel &Model, size_t freqindex,
    std::string TempDirName);

HPX_DEFINE_PLAIN_ACTION(HPXCalculateFrequency, CalculateFrequency_action);


#endif /* CALCFREQ_H_ */
