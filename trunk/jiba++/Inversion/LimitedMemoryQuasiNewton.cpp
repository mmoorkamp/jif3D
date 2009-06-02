//============================================================================
// Name        : LimitedMemoryQuasiNewton.cpp
// Author      : Apr 20, 2009
// Version     :
// Copyright   : 2009, mmoorkamp
//============================================================================

#include "../Global/FatalException.h"
#include "LimitedMemoryQuasiNewton.h"
#include "BacktrackingLineSearch.h"
#include "mcsrch.h"
#include <fstream>
namespace jiba
  {

    LimitedMemoryQuasiNewton::LimitedMemoryQuasiNewton(boost::shared_ptr<
        jiba::ObjectiveFunction> ObjFunction, const size_t n) :
      NonLinearOptimization(ObjFunction), mu(1),MaxPairs(n)
      {

      }

    LimitedMemoryQuasiNewton::~LimitedMemoryQuasiNewton()
      {

      }

    void LimitedMemoryQuasiNewton::StepImplementation(jiba::rvec &CurrentModel)
      {

        SearchDir = CovGrad;
        const size_t nmod = CovGrad.size();
        const size_t npairs = SHistory.size();

        jiba::rvec Alpha(npairs), Rho(npairs);

        //we store the elements in reverse order
        for (int i = npairs - 1; i >= 0; --i)
          {
            Rho( i) = 1. / ublas::inner_prod(*YHistory.at(i),
                ublas::element_div(*SHistory.at(i), GetModelCovDiag()));
            Alpha( i) = Rho(i) * ublas::inner_prod(*SHistory.at(i),
                ublas::element_div(SearchDir, GetModelCovDiag()));
            SearchDir -= Alpha(i) * *YHistory.at(i);
          }

        for (size_t i = 0; i < npairs; ++i)
          {
            double beta = Rho(i) * ublas::inner_prod(*YHistory.at(i),
                ublas::element_div(SearchDir, GetModelCovDiag()));
            SearchDir += *SHistory.at(i) * (Alpha(i) - beta);
          }
        SearchDir *= -1.0;
        //std::cout << "Raw Gradient: " << RawGrad << std::endl;
        //std::cout << "Gradient: " << Gradient << std::endl;
        //std::cout << "Search Dir: " << SearchDir << std::endl;
        //double mu = BacktrackingLineSearch().FindStep(CurrentModel, Gradient,
        //    SearchDir, *GetObjective());
        //std::cout << "SearchDir: " << SearchDir << std::endl;
        //std::cout << "RawGrad: " << RawGrad << std::endl;

        int status = OPTPP::mcsrch(GetObjective().get(), SearchDir, RawGrad, CurrentModel, Misfit,
            &mu, 20, 1e-4, 2.2e-16, 0.9, 1e9, 1e-9);
        if (status < 0)
        	throw jiba::FatalException("Cannot find suitable step");
        std::cout << " Mu: " << mu  << std::endl;
        CurrentModel += mu * SearchDir;
        if (npairs < MaxPairs)
          {
            boost::shared_ptr<jiba::rvec> NewS(new jiba::rvec(nmod));
            SHistory.push_back(NewS);
            boost::shared_ptr<jiba::rvec> NewY(new jiba::rvec(nmod));
            YHistory.push_back(NewY);
          }
        else
          {
            std::rotate(SHistory.begin(), SHistory.begin() + 1, SHistory.end());
            std::rotate(YHistory.begin(), YHistory.begin() + 1, YHistory.end());
          }
        *SHistory.back() = mu * SearchDir;
        *YHistory.back() = ublas::element_prod(RawGrad, GetModelCovDiag()) - CovGrad;
        CovGrad = ublas::element_prod(RawGrad, GetModelCovDiag());
      }
  }
