//============================================================================
// Name        : BacktrackingLineSearch.h
// Author      : Apr 20, 2009
// Version     :
// Copyright   : 2009, mmoorkamp
//============================================================================


#ifndef BACKTRACKINGLINESEARCH_H_
#define BACKTRACKINGLINESEARCH_H_

#include "LineSearch.h"

namespace jiba
  {

    class BacktrackingLineSearch: public jiba::LineSearch
      {
    public:
      virtual double FindStep(const jiba::rvec &CurrModel,const jiba::rvec &CurrGradient,const jiba::rvec &CurrSearch,ObjectiveFunction &Objective);
      BacktrackingLineSearch();
      virtual ~BacktrackingLineSearch();
      };

  }

#endif /* BACKTRACKINGLINESEARCH_H_ */
