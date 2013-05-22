//============================================================================
// Name        : NonLinearConjugateGradient.cpp
// Author      : Apr 20, 2009
// Version     :
// Copyright   : 2009, mmoorkamp
//============================================================================

#include "../Global/FatalException.h"
#include "../Global/convert.h"
#include "../Global/NormProd.h"
#include "NonLinearConjugateGradient.h"
#include "mcsrch.h"

namespace jif3D
  {

    NonLinearConjugateGradient::NonLinearConjugateGradient(boost::shared_ptr<
        jif3D::ObjectiveFunction> ObjFunction) :
      GradientBasedOptimization(ObjFunction), OldGradient(), OldDirection(),
          OldOmega(1.0), mu(1.0)
      {

      }

    NonLinearConjugateGradient::~NonLinearConjugateGradient()
      {

      }

    void NonLinearConjugateGradient::StepImplementation(
        jif3D::rvec &CurrentModel)
      {
        const size_t nmod = CovGrad.size();
        //if we are in the first iteration or we took
        //a steepest descent step, reset all information
        if (OldGradient.size() != nmod)
          {
            OldOmega = 1.0;
            OldGradient.resize(nmod);
            OldDirection.resize(nmod);
            OldGradient.clear();
            OldDirection.clear();
            gamma = 1.0;
          }
        double omega = 0.0;
        double alpha = 0.0;
        //doing it in a loop instead of using inner_product and element_div
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
        SearchDir *= -gamma;
        //ublas::element_div(SearchDir, M);
        //we only care about the sign of the angle, so we don't bother normalizing
        //by the norm of the two vectors
        double angle = ublas::inner_prod(SearchDir, RawGrad);
        int status = 0;
        //if the search direction is a descent direction use it
        //otherwise use the direction of steepest descent
        if (angle < 0.0)
          {
            status = OPTPP::mcsrch(&GetObjective(), SearchDir, RawGrad,
                CurrentModel, Misfit, &mu, 20, 1e-4, 2.2e-16, 0.1, 1e9, 1e-9);
            jif3D::rvec y(ublas::element_prod(RawGrad, GetModelCovDiag())
                - OldGradient);
            //This is the same calculation as for L-BFGS to scale the search direction for the next iteration
            //we therefore use the same nomenclature
            jif3D::rvec s(mu * SearchDir);
            double rho = 1.0 / NormProd(y, s, GetModelCovDiag());
            gamma = 1.0 / rho / NormProd(y, y, GetModelCovDiag());
          }
        else
          {
            CovGrad *= -1.0;
            OldGradient.resize(0);
            status = OPTPP::mcsrch(&GetObjective(), CovGrad, RawGrad,
                CurrentModel, Misfit, &mu, 20, 1e-4, 2.2e-16, 0.1, 1e9, 1e-9);
          }
        CurrentModel += mu * SearchDir;
        //the line search already calculated the new gradient
        // so we signal that we have calculated everything
        HaveEvaluated();
      }
  }
