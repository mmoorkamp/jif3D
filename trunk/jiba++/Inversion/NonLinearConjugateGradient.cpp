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

        std::cout << "Misfit: " << Misfit << " Model: " << CurrentModel
            << " Gradient: " << ublas::norm_2(Gradient) << std::endl;
        const size_t nmod = Gradient.size();
        if (OldGradient.size() != nmod)
          {
            OldOmega = 1.0;
            OldGradient.resize(nmod);
            OldDirection.resize(nmod);
            std::fill(OldGradient.begin(), OldGradient.end(), 0.0);
            std::fill(OldDirection.begin(), OldDirection.end(), 0.0);
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
        jiba::rvec Search(Gradient + alpha * OldDirection);
        OldDirection = Search;

        Search *= -1.0;
        mu = 1.0;
        double angle = ublas::inner_prod(Search, RawGrad);
        std::cout << "Angle: " << angle << std::endl;
        if (angle < 0.0)
          {
            OPTPP::mcsrch(GetObjective().get(), Search, RawGrad, CurrentModel,
                Misfit, &mu, 20, 1e-4, 2.2e-16, 0.1, 1e9, 1e-9);
          }
        else
          {
            Gradient *= -1.0;
            OldGradient.resize(0);
            OPTPP::mcsrch(GetObjective().get(), Gradient, RawGrad,
                CurrentModel, Misfit, &mu, 20, 1e-4, 2.2e-16, 0.1, 1e9, 1e-9);
          }

        CurrentModel += mu * Search;

      }
  }
