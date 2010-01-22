//============================================================================
// Name        : LimitedMemoryQuasiNewton.cpp
// Author      : Apr 20, 2009
// Version     :
// Copyright   : 2009, mmoorkamp
//============================================================================

#include "../Global/FatalException.h"
#include "../Global/convert.h"
#include "LimitedMemoryQuasiNewton.h"
#include "mcsrch.h"
namespace jiba
  {

    LimitedMemoryQuasiNewton::LimitedMemoryQuasiNewton(boost::shared_ptr<
        jiba::ObjectiveFunction> ObjFunction, const size_t n) :
      GradientBasedOptimization(ObjFunction), mu(1.0), LineIter(20),MaxPairs(n), SHistory(),
          YHistory()
      {

      }

    LimitedMemoryQuasiNewton::~LimitedMemoryQuasiNewton()
      {

      }

    void LimitedMemoryQuasiNewton::StepImplementation(jiba::rvec &CurrentModel)
      {
        //we start assuming the search direction is the direction of steepest ascent
        SearchDir = CovGrad;
        const size_t nmod = CovGrad.size();
        const size_t npairs = SHistory.size();

        jiba::rvec Alpha(npairs), Rho(npairs);

        //we store the elements in reverse order
        //and apply algorithm 9.1 from Nocedal and Wright
        for (int i = npairs - 1; i >= 0; --i)
          {
            Rho(i) = 1. / ublas::inner_prod(*YHistory.at(i),
                ublas::element_div(*SHistory.at(i), GetModelCovDiag()));
            Alpha(i) = Rho(i) * ublas::inner_prod(*SHistory.at(i),
                ublas::element_div(SearchDir, GetModelCovDiag()));
            SearchDir -= Alpha(i) * *YHistory.at(i);
          }
        double gamma = 1.0;
        if (YHistory.size() > 0)
          {
            gamma = 1.0 / Rho(npairs-1) / ublas::inner_prod(*YHistory.back(),
                ublas::element_div(*YHistory.back(), GetModelCovDiag()));
          }
        SearchDir *= gamma;
        for (size_t i = 0; i < npairs; ++i)
          {
            double beta = Rho(i) * ublas::inner_prod(*YHistory.at(i),
                ublas::element_div(SearchDir, GetModelCovDiag()));
            SearchDir += *SHistory.at(i) * (Alpha(i) - beta);
          }
        //at each iteration we reset the stepsize
        mu = 1.0;
        SearchDir *= -1.0;
        //now we do a line search to find the optimum step size mu
        //after this call, both Misfit and RawGrad are already
        //updated for the new model
        int status = OPTPP::mcsrch(&GetObjective(), SearchDir, RawGrad,
            CurrentModel, Misfit, &mu, LineIter, 1e-4, 2.2e-16, 0.9, 1e9, 1e-12);

        if (status < 0 )
          {
            throw jiba::FatalException("Cannot find suitable step. Status: "
                + jiba::stringify(status));
          }
        //if we have found a good stepsize, update the model
        CurrentModel += mu * SearchDir;
        //if we haven't reached the maximum number of correction pairs, yet
        if (npairs < MaxPairs)
          {
            //allocate storage for a new correction pair
            boost::shared_ptr<jiba::rvec> NewS(new jiba::rvec(nmod));
            SHistory.push_back(NewS);
            boost::shared_ptr<jiba::rvec> NewY(new jiba::rvec(nmod));
            YHistory.push_back(NewY);
          }
        else
          {
            //otherwise shift the correction pairs in the storage
            //so the first (oldest) becomes the last element
            std::rotate(SHistory.begin(), SHistory.begin() + 1, SHistory.end());
            std::rotate(YHistory.begin(), YHistory.begin() + 1, YHistory.end());
          }
        //and overwrite the last (oldest) correction pair
        *SHistory.back() = mu * SearchDir;
        *YHistory.back() = ublas::element_prod(RawGrad, GetModelCovDiag())
            - CovGrad;
        //the line search has updated the gradient and misfit
        HaveEvaluated();
      }
  }
