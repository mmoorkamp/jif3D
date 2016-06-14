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
        jif3D::ObjectiveFunction> ObjFunction) :
      ModelCovDiag(), Objective(ObjFunction), Misfit()
      {

      }

    NonLinearOptimization::~NonLinearOptimization()
      {

      }

    void NonLinearOptimization::MakeStep(jif3D::rvec &CurrentModel)
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
