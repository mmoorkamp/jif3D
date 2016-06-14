//============================================================================
// Name        : ModEMCalculator.cpp
// Author      : 14 Jan 2016
// Version     : 
// Copyright   : 2016, mm489
//============================================================================

#include "ModEMCalculator.h"

namespace jif3D
  {

    ModEMCalculator::ModEMCalculator()
      {
      }

    ModEMCalculator::~ModEMCalculator()
      {
      }

    rvec ModEMCalculator::Calculate(const ModelType &Model, size_t minfreqindex,
        size_t maxfreqindex)
      {
    	return rvec();
      }

    rvec ModEMCalculator::LQDerivative(const ModelType &Model, const rvec &Misfit,
        size_t minfreqindex, size_t maxfreqindex)
      {
    	return rvec();
      }

  } /* namespace jif3D */
