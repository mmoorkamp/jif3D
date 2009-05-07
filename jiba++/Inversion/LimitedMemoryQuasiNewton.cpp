//============================================================================
// Name        : LimitedMemoryQuasiNewton.cpp
// Author      : Apr 20, 2009
// Version     :
// Copyright   : 2009, mmoorkamp
//============================================================================


#include "LimitedMemoryQuasiNewton.h"
#include "BacktrackingLineSearch.h"
#include "mcsrch.h"
#include <fstream>
namespace jiba
  {

    LimitedMemoryQuasiNewton::LimitedMemoryQuasiNewton(boost::shared_ptr<
        jiba::ObjectiveFunction> ObjFunction, const size_t n) :
      NonLinearOptimization(ObjFunction), MaxPairs(n)
      {

      }

    LimitedMemoryQuasiNewton::~LimitedMemoryQuasiNewton()
      {

      }

    void LimitedMemoryQuasiNewton::StepImplementation(jiba::rvec &CurrentModel)
      {
        Misfit = GetObjective()->CalcMisfit(CurrentModel);
        std::cout << "Misfit: " << Misfit << std::endl;
        jiba::rvec RawGrad(GetObjective()->CalcGradient());
        jiba::rvec Gradient(ublas::element_prod(RawGrad, GetModelCovDiag()));

        const size_t nmod = Gradient.size();
        const size_t npairs = SHistory.size();

        jiba::rvec Alpha(npairs), Rho(npairs), SearchDir(Gradient);

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

        //double mu = BacktrackingLineSearch().FindStep(CurrentModel, Gradient,
        //    SearchDir, *GetObjective());
        //std::cout << "SearchDir: " << SearchDir << std::endl;
        //std::cout << "RawGrad: " << RawGrad << std::endl;
        double mu = 1.0;
        int status = OPTPP::mcsrch(GetObjective().get(), SearchDir, RawGrad, CurrentModel, Misfit,
            &mu, 20, 1e-4, 2.2e-16, 0.9, 1e9, 1e-9);
        std::cout << "Status: " << status << std::endl;
        std::cout << " Mu: " << mu << "Search Angle: " << ublas::inner_prod(
                   SearchDir, RawGrad) / (ublas::norm_2(SearchDir) * ublas::norm_2(
                   RawGrad)) << std::endl;
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
        *YHistory.back() = ublas::element_prod(GetObjective()->CalcGradient(), GetModelCovDiag()) - Gradient;
      }
  }
