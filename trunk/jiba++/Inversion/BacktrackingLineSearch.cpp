//============================================================================
// Name        : BacktrackingLineSearch.cpp
// Author      : Apr 20, 2009
// Version     :
// Copyright   : 2009, mmoorkamp
//============================================================================


#include "BacktrackingLineSearch.h"
#include "../Global/VecMat.h"
namespace jiba
  {

    BacktrackingLineSearch::BacktrackingLineSearch()
      {

      }

    BacktrackingLineSearch::~BacktrackingLineSearch()
      {

      }

    double BacktrackingLineSearch::FindStep(const jiba::rvec &CurrModel,
        const jiba::rvec &CurrGradient, const jiba::rvec &CurrSearch,
        ObjectiveFunction &Objective)
      {
        double rho = 0.5;
        double c = 1e-5;
        double alpha = 1.0;

        const double Misfit = Objective.CalcMisfit(CurrModel);
        jiba::rvec TrialModel = CurrModel + alpha * CurrSearch;
        double DecFactor = c * boost::numeric::ublas::inner_prod(CurrGradient,CurrSearch);
        //std::cout << "DecFactor: " << DecFactor << std::endl;
        double NewMisfit = Objective.CalcMisfit(TrialModel);
        while (NewMisfit > Misfit - alpha * DecFactor && alpha > 1e-10)
          {
            alpha *= rho;
            TrialModel = CurrModel + alpha * CurrSearch;
            NewMisfit = Objective.CalcMisfit(TrialModel);
            //std::cout << "Alpha: " << alpha << " Trial Model: " << TrialModel << " NewMisfit: " << NewMisfit << std::endl;
          }

        return alpha;
      }
  }
