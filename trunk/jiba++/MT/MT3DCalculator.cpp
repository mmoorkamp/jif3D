//============================================================================
// Name        : MT3DCalculator.cpp
// Author      : Jul 6, 2009
// Version     : 
// Copyright   : 2009, mmoorkamp
//============================================================================


#include "MT3DCalculator.h"

namespace jiba
  {

    rvec MT3DCalculator::Calculate(const ThreeDMTModel &Model)
      {
        return CalculateImpl(Model);
      }
    rvec MT3DCalculator::LQDerivative(const ThreeDMTModel &Model,
        const rvec &Misfit)
      {
        return LQDerivativeImpl(Model, Misfit);
      }

    MT3DCalculator::MT3DCalculator()
      {

      }

    MT3DCalculator::~MT3DCalculator()
      {

      }

  }
