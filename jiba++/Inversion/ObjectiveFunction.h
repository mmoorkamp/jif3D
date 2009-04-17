//============================================================================
// Name        : ObjectiveFunction.h
// Author      : Apr 16, 2009
// Version     :
// Copyright   : 2009, mmoorkamp
//============================================================================


#ifndef OBJECTIVEFUNCTION_H_
#define OBJECTIVEFUNCTION_H_

#include "../Global/VecMat.h"

namespace jiba
  {
    /** \addtogroup inversion General routines for inversion */
    /* @{ */
    class ObjectiveFunction
      {
    private:
      jiba::rvec DataDifference;
      virtual void
          ImplDataDifference(const jiba::rvec &Model, jiba::rvec &Diff) = 0;
      virtual jiba::rvec ImplGradient(const jiba::rvec &Model, const  jiba::rvec &Diff) = 0;
    public:
      double CalcMisfit(const jiba::rvec &Model)
        {
          ImplDataDifference(Model, DataDifference);
          return boost::numeric::ublas::norm_2(DataDifference);
        }
      jiba::rvec CalcGradient(const jiba::rvec &Model)
        {
          return ImplGradient(Model, DataDifference);
        }
      ObjectiveFunction();
      virtual ~ObjectiveFunction();
      };
    /* @} */
  }

#endif /* OBJECTIVEFUNCTION_H_ */
