//============================================================================
// Name        : ObjectiveFunction.h
// Author      : Apr 16, 2009
// Version     :
// Copyright   : 2009, mmoorkamp
//============================================================================


#ifndef OBJECTIVEFUNCTION_H_
#define OBJECTIVEFUNCTION_H_

#include "../Global/VecMat.h"
#include <cassert>

namespace jiba
  {
    /** \addtogroup inversion General routines for inversion */
    /* @{ */
    class ObjectiveFunction
      {
    private:
      jiba::rvec DataDifference;
      jiba::rvec CovarDiag;
      virtual void
          ImplDataDifference(const jiba::rvec &Model, jiba::rvec &Diff) = 0;
      virtual jiba::rvec ImplGradient(const jiba::rvec &Model,
          const jiba::rvec &Diff) = 0;
    public:
      void SetDataCovar(const jiba::rvec &Cov)
        {
          CovarDiag = Cov;
        }
      double CalcMisfit(const jiba::rvec &Model)
        {
          ImplDataDifference(Model, DataDifference);
          if (CovarDiag.size() != DataDifference.size())
            {
              CovarDiag.resize(DataDifference.size());
              std::fill(CovarDiag.begin(), CovarDiag.end(), 1.0);
            }
          return ublas::inner_prod(DataDifference, ublas::element_div(
              DataDifference, CovarDiag));
        }
      jiba::rvec CalcGradient(const jiba::rvec &Model)
        {
          assert(CovarDiag.size() == DataDifference.size());
          return ImplGradient(Model, ublas::element_div(DataDifference,
              CovarDiag));
        }
      ObjectiveFunction();
      virtual ~ObjectiveFunction();
      };
  /* @} */
  }

#endif /* OBJECTIVEFUNCTION_H_ */
