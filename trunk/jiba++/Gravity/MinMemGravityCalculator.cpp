//============================================================================
// Name        : MinMemGravityCalculator.cpp
// Author      : Dec 10, 2008
// Version     :
// Copyright   : 2008, mmoorkamp
//============================================================================


#include "MinMemGravityCalculator.h"

namespace jiba
  {

    MinMemGravityCalculator::MinMemGravityCalculator(ThreeDGravityImplementation &TheImp):
      ThreeDGravityCalculator(TheImp)
      {

      }

    MinMemGravityCalculator::~MinMemGravityCalculator()
      {
        // TODO Auto-generated destructor stub
      }

    rvec MinMemGravityCalculator::Calculate(const ThreeDGravityModel &Model)
      {
    	SetCurrentSensitivities().resize(0,0);
        CheckModelConsistency(Model);

        return Imp.Calculate(Model,*this);

      }

  }
