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
      ModelCovDiag(), Objective(ObjFunction), LastModel(), RawGrad(),
          CovGrad(), SearchDir(), Misfit()
      {

      }

    NonLinearOptimization::~NonLinearOptimization()
      {

      }

    void NonLinearOptimization::CalcMisfitandGradient(
        const jiba::rvec &CurrentModel)
      {
        Misfit = GetObjective()->CalcMisfit(CurrentModel);

        RawGrad = GetObjective()->CalcGradient();
        CovGrad = ublas::element_prod(RawGrad, GetModelCovDiag());
      }

    void NonLinearOptimization::MakeStep(jiba::rvec &CurrentModel)
      {
        if (ModelCovDiag.size() != CurrentModel.size())
          {
            ModelCovDiag.resize(CurrentModel.size());
            std::fill(ModelCovDiag.begin(), ModelCovDiag.end(), 1.0);
          }
        if (SearchDir.size() != CurrentModel.size())
          {
            SearchDir.resize(CurrentModel.size());
          }
        if (LastModel.size() != CurrentModel.size() || !std::equal(
            LastModel.begin(), LastModel.end(), CurrentModel.begin()))
          {
            CalcMisfitandGradient(CurrentModel);
          }
        StepImplementation(CurrentModel);
        LastModel = CurrentModel;
      }
  }
