//============================================================================
// Name        : NonLinearConjugateGradient.cpp
// Author      : Apr 20, 2009
// Version     :
// Copyright   : 2009, mmoorkamp
//============================================================================


#include "NonLinearConjugateGradient.h"
#include "BacktrackingLineSearch.h"
#include "mcsrch.h"
#include <fstream>
namespace jiba
  {

    NonLinearConjugateGradient::NonLinearConjugateGradient(boost::shared_ptr<
        jiba::ObjectiveFunction> ObjFunction) :
      NonLinearOptimization(ObjFunction), OldOmega(1.0), mu(1.0)
      {

      }

    NonLinearConjugateGradient::~NonLinearConjugateGradient()
      {

      }

    void NonLinearConjugateGradient::StepImplementation(
        jiba::rvec &CurrentModel)
      {
        Misfit = GetObjective()->CalcMisfit(CurrentModel);

        jiba::rvec RawGrad(GetObjective()->CalcGradient());
        jiba::rvec Gradient(ublas::element_prod(RawGrad, GetModelCovDiag()));
        std::cout << "Misfit: " << Misfit << " Gradient: " << ublas::norm_2(
            Gradient) << std::endl;
        const size_t nmod = Gradient.size();
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
        for (size_t i = 0; i < nmod; ++i)
          {
            const double factor = 1.0 / GetModelCovDiag()(i) * Gradient(i);
            omega += Gradient(i) * factor;
            alpha += OldGradient(i) * factor;
          }
        alpha = (omega - alpha) / OldOmega;
        OldGradient = Gradient;
        OldOmega = omega;
        SearchDir = Gradient + alpha * OldDirection;
        OldDirection = SearchDir;

        SearchDir *= -1.0;
        mu = 1.0;
        double angle = ublas::inner_prod(SearchDir, RawGrad) / (ublas::norm_2(
            SearchDir) * ublas::norm_2(RawGrad));
        std::cout << "Angle: " << angle << std::endl;
        int status = 0;
        if (angle < 0.0)
          {
            status = OPTPP::mcsrch(GetObjective().get(), SearchDir, RawGrad,
                CurrentModel, Misfit, &mu, 20, 1e-4, 2.2e-16, 0.1, 1e9, 1e-9);
          }
        else
          {
            Gradient *= -1.0;
            OldGradient.resize(0);
            status = OPTPP::mcsrch(GetObjective().get(), Gradient, RawGrad,
                CurrentModel, Misfit, &mu, 20, 1e-4, 2.2e-16, 0.1, 1e9, 1e-9);
          }
        std::cout << "Status: " << status << std::endl;
        CurrentModel += mu * SearchDir;

      }
  }
