//============================================================================
// Name        : LimitedMemoryQuasiNewton.cpp
// Author      : Apr 20, 2009
// Version     :
// Copyright   : 2009, mmoorkamp
//============================================================================

#include "../Global/FatalException.h"
#include "../Global/convert.h"
#include "../Global/NormProd.h"
#include "LimitedMemoryQuasiNewton.h"
#include "mcsrch.h"
#include <iostream>

namespace jif3D
  {

    LimitedMemoryQuasiNewton::LimitedMemoryQuasiNewton(
        boost::shared_ptr<jif3D::ObjectiveFunction> ObjFunction, boost::shared_ptr<jif3D::GeneralCovariance> Cv, const size_t n,
        bool gradscale) :
        GradientBasedOptimization(ObjFunction, Cv), scale(gradscale), mu(1.0), LineIter(10), MaxPairs(
            n), SHistory(), YHistory()
      {

      }

    LimitedMemoryQuasiNewton::~LimitedMemoryQuasiNewton()
      {

      }

    void LimitedMemoryQuasiNewton::StepImplementation(jif3D::rvec &CurrentModel)
      {
        //we start assuming the search direction is the direction of steepest ascent
        SearchDir = CovGrad;
        const size_t nmod = CovGrad.size();
        const size_t npairs = SHistory.size();

        jif3D::rvec Alpha(npairs), Rho(npairs);

        //we store the elements in reverse order
        //and apply algorithm 9.1 from Nocedal and Wright
        for (int i = npairs - 1; i >= 0; --i)
          {
            Rho(i) = 1. / NormProd(*YHistory.at(i), *SHistory.at(i), GetCovObj());
            Alpha(i) = Rho(i) * NormProd(*SHistory.at(i), SearchDir, GetCovObj());
            SearchDir -= Alpha(i) * *YHistory.at(i);
          }
        //gamma is a scaling factor that we apply to the search direction
        //after the first iteration this scale usually helps to find the minimum
        //without additional steps in the line search
        double gamma = 1.0;
        if (npairs > 0)
          {
            gamma = 1.0 / Rho(npairs - 1)
                / NormProd(*YHistory.back(), *YHistory.back(), GetCovObj());
          }
        else
          {
            //if we don't have the previous history, we use the length of the
            //gradient as a scale
            if (scale)
              {
                //we only get here the first time when searchdir is equal to CovGrad
                // and algorithm 9.1 above has not worked yet
                //theoretically the norm should use the search direction and the covariance as norm
                //however, we can save the multiplication by the model covariance and its inverse
                gamma = 1.0 / sqrt(ublas::inner_prod(RawGrad, RawGrad));
              }
          }
        SearchDir *= -gamma;
        for (size_t i = 0; i < npairs; ++i)
          {
            double beta = Rho(i)
                * NormProd(*YHistory.at(i), SearchDir, GetCovObj());
            SearchDir -= *SHistory.at(i) * (Alpha(i) + beta);
          }
        //at each iteration we reset the stepsize
        mu = 1.0;
        //now we do a line search to find the optimum step size mu
        //after this call, both Misfit and RawGrad are already
        //updated for the new model
        int status = OPTPP::mcsrch(&GetObjective(), SearchDir, RawGrad, CurrentModel,
            Misfit, &mu, LineIter, 1e-4, 2.2e-16, 0.9, 1e9, 1e-12);

        if (status < 0)
          {
            std::cerr << "Cannot find suitable step. Status: " + jif3D::stringify(status)
                << std::endl;
            throw jif3D::FatalException(
                "Cannot find suitable step. Status: " + jif3D::stringify(status),
                __FILE__, __LINE__);
          }
        //if we have found a good stepsize, update the model
        CurrentModel += mu * SearchDir;
        //if we haven't reached the maximum number of correction pairs, yet
        if (npairs < MaxPairs)
          {
            //allocate storage for a new correction pair
            boost::shared_ptr<jif3D::rvec> NewS(new jif3D::rvec(nmod));
            SHistory.push_back(NewS);
            boost::shared_ptr<jif3D::rvec> NewY(new jif3D::rvec(nmod));
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
        *YHistory.back() = GetCovObj()->ApplyCovar(RawGrad) - CovGrad;
        //the line search has updated the gradient and misfit
        HaveEvaluated();
      }
  }
