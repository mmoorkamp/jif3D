//============================================================================
// Name        : ObjectiveFunction.cpp
// Author      : Apr 16, 2009
// Version     :
// Copyright   : 2009, mmoorkamp
//============================================================================

#include "ObjectiveFunction.h"
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>

namespace jif3D
  {

    ObjectiveFunction::ObjectiveFunction() :
        nEval(0), DataDifference(), InvCovMat()
      {

      }

    ObjectiveFunction::~ObjectiveFunction()
      {

      }

  }
//BOOST_CLASS_EXPORT_IMPLEMENT(jif3D::ObjectiveFunction)
