//============================================================================
// Name        : NonLinearOptimization.cpp
// Author      : Apr 16, 2009
// Version     :
// Copyright   : 2009, mmoorkamp
//============================================================================


#include "NonLinearOptimization.h"

namespace jiba
  {

    NonLinearOptimization::NonLinearOptimization(boost::shared_ptr<
        jiba::ObjectiveFunction> ObjFunction) :
      ModelCovDiag(), Objective(ObjFunction), Misfit()
      {

      }

    NonLinearOptimization::~NonLinearOptimization()
      {

      }

    void NonLinearOptimization::MakeStep(jiba::rvec &CurrentModel)
      {
        if (ModelCovDiag.size() != CurrentModel.size())
          {
            ModelCovDiag.resize(CurrentModel.size());
            std::fill(ModelCovDiag.begin(), ModelCovDiag.end(), 1.0);
          }

        EvaluateModel(CurrentModel);
        StepImplementation(CurrentModel);
      }
  }
