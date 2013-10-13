//============================================================================
// Name        : MinMemGravityCalculator.cpp
// Author      : Dec 10, 2008
// Version     :
// Copyright   : 2008, mmoorkamp
//============================================================================


#include "MinMemGravMagCalculator.h"
#include <iostream>


namespace jif3D
  {

    MinMemGravMagCalculator::MinMemGravMagCalculator(boost::shared_ptr<
        ThreeDGravMagImplementation> TheImp) :
      ThreeDGravMagCalculator(TheImp)
      {

      }

    MinMemGravMagCalculator::~MinMemGravMagCalculator()
      {

      }

    /*! For this simple case the calculator class only checks the model for consistency,
     * then calls the implementation object to do the calculation and returns the result
     * @param Model The Gravity model for the forward calculation
     * @return The resulting measurements (scalar or tensorial)
     */
    rvec MinMemGravMagCalculator::Calculate(const ThreeDGravityModel &Model)
      {
        SetCurrentSensitivities().resize(0, 0);
        return ThreeDGravMagCalculator::Calculate(Model);
      }

  }
