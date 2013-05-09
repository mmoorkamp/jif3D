//============================================================================
// Name        : MinMemGravityCalculator.cpp
// Author      : Dec 10, 2008
// Version     :
// Copyright   : 2008, mmoorkamp
//============================================================================


#include "MinMemGravityCalculator.h"
#include <iostream>


namespace jiba
  {

    MinMemGravityCalculator::MinMemGravityCalculator(boost::shared_ptr<
        ThreeDGravityImplementation> TheImp) :
      ThreeDGravityCalculator(TheImp)
      {

      }

    MinMemGravityCalculator::~MinMemGravityCalculator()
      {

      }

    /*! For this simple case the calculator class only checks the model for consistency,
     * then calls the implementation object to do the calculation and returns the result
     * @param Model The Gravity model for the forward calculation
     * @return The resulting measurements (scalar or tensorial)
     */
    rvec MinMemGravityCalculator::Calculate(const ThreeDGravityModel &Model)
      {
        SetCurrentSensitivities().resize(0, 0);
        return ThreeDGravityCalculator::Calculate(Model);
      }

  }
