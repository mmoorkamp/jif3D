//============================================================================
// Name        : LineSearch.h
// Author      : Apr 20, 2009
// Version     :
// Copyright   : 2009, mmoorkamp
//============================================================================


#ifndef LINESEARCH_H_
#define LINESEARCH_H_

#include "../Global/VecMat.h"
#include "ObjectiveFunction.h"

namespace jiba
  {
    /** \addtogroup inversion General routines for inversion */
    /* @{ */
    class LineSearch
      {
    public:
      virtual double FindStep(const jiba::rvec &CurrModel,
          const jiba::rvec &CurrGradient, const jiba::rvec &CurrSearch,
          ObjectiveFunction &Objective) = 0;
      LineSearch();
      virtual ~LineSearch();
      };
  /* @} */
  }

#endif /* LINESEARCH_H_ */
