//============================================================================
// Name        : ObjectiveFunction.cpp
// Author      : Apr 16, 2009
// Version     :
// Copyright   : 2009, mmoorkamp
//============================================================================


#include "ObjectiveFunction.h"

namespace jiba
  {

    ObjectiveFunction::ObjectiveFunction() :
      nEval(0), DataDifference(), CovarDiag(), PreCondDiag(),
          DataTransform()
      {

      }

    ObjectiveFunction::~ObjectiveFunction()
      {

      }

  }