//============================================================================
// Name        : mcsrch.h
// Author      : Apr 20, 2009
// Version     :
// Copyright   : 2009, mmoorkamp
//============================================================================


#ifndef MCSRCH_H_
#define MCSRCH_H_
#include "ObjectiveFunction.h"
#include "../Global/VecMat.h"
namespace OPTPP
  {
    int mcsrch(jiba::ObjectiveFunction* nlp, const jiba::rvec& s,
        const jiba::rvec &model, double misfit, double *stp, int itnmax, double ftol,
        double xtol, double gtol, double stpmax, double stpmin);
  }
#endif /* MCSRCH_H_ */
