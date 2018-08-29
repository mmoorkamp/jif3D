//============================================================================
// Name        : NonLinearOptimization.cpp
// Author      : Apr 16, 2009
// Version     :
// Copyright   : 2009, mmoorkamp
//============================================================================


#include "NonLinearOptimization.h"

namespace jif3D
  {

    NonLinearOptimization::NonLinearOptimization(boost::shared_ptr<
        jif3D::ObjectiveFunction> ObjFunction, boost::shared_ptr<jif3D::GeneralCovariance> Cv) :
      Covar(Cv), Objective(ObjFunction), Misfit()
      {

      }

    NonLinearOptimization::~NonLinearOptimization()
      {

      }

    void NonLinearOptimization::MakeStep(jif3D::rvec &CurrentModel)
      {
        EvaluateModel(CurrentModel);
        StepImplementation(CurrentModel);
      }
  }
