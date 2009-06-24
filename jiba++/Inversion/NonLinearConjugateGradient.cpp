//============================================================================
// Name        : NonLinearConjugateGradient.cpp
// Author      : Apr 20, 2009
// Version     :
// Copyright   : 2009, mmoorkamp
//============================================================================


#include "NonLinearConjugateGradient.h"
#include "mcsrch.h"

namespace jiba
  {

    NonLinearConjugateGradient::NonLinearConjugateGradient(boost::shared_ptr<
        jiba::ObjectiveFunction> ObjFunction) :
      NonLinearOptimization(ObjFunction), OldGradient(),
          OldDirection(), OldOmega(1.0), mu(1.0)
      {

      }

    NonLinearConjugateGradient::~NonLinearConjugateGradient()
      {

      }

    void NonLinearConjugateGradient::StepImplementation(
        jiba::rvec &CurrentModel)
      {
        const size_t nmod = CovGrad.size();
        if (OldGradient.size() != nmod)
          {
            OldOmega = 1.0;
            OldGradient.resize(nmod);
            OldDirection.resize(nmod);
            OldGradient.clear();
            OldDirection.clear();
          }
        double omega = 0.0;
        double alpha = 0.0;
        //doing it like this instead of using inner_product and element_div
        //saves us a temporary object and is faster
        for (size_t i = 0; i < nmod; ++i)
          {
            const double factor = CovGrad(i) / GetModelCovDiag()(i);
            omega += CovGrad(i) * factor;
            alpha += OldGradient(i) * factor;
          }
        alpha = (omega - alpha) / OldOmega;
        OldGradient = CovGrad;
        OldOmega = omega;
        SearchDir = CovGrad + alpha * OldDirection;
        OldDirection = SearchDir;
        mu = 1.0;
        SearchDir *= -1.0;
        double angle = ublas::inner_prod(SearchDir, RawGrad) / (ublas::norm_2(
            SearchDir) * ublas::norm_2(RawGrad));
        int status = 0;
        if (angle < 0.0)
          {
            status = OPTPP::mcsrch(GetObjective().get(), SearchDir, RawGrad,
                CurrentModel, Misfit, &mu, 20, 1e-4, 2.2e-16, 0.1, 1e9, 1e-9);
          }
        else
          {
            CovGrad *= -1.0;
            OldGradient.resize(0);
            status = OPTPP::mcsrch(GetObjective().get(), CovGrad, RawGrad,
                CurrentModel, Misfit, &mu, 20, 1e-4, 2.2e-16, 0.1, 1e9, 1e-9);
          }
        CurrentModel += mu * SearchDir;
        CovGrad = ublas::element_prod(RawGrad, GetModelCovDiag());
      }
  }
