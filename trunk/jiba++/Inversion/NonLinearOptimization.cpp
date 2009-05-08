//============================================================================
// Name        : NonLinearOptimization.cpp
// Author      : Apr 16, 2009
// Version     :
// Copyright   : 2009, mmoorkamp
//============================================================================


#include "NonLinearOptimization.h"

namespace jiba
  {

    NonLinearOptimization::NonLinearOptimization(boost::shared_ptr<jiba::ObjectiveFunction> ObjFunction)
    :Objective(ObjFunction)
      {

      }

    NonLinearOptimization::~NonLinearOptimization()
      {

      }

  }
