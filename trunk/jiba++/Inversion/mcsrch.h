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
    //! A version of the More and Thuente line search algorithm taken from Opt++
    int mcsrch(jiba::ObjectiveFunction* nlp, const jiba::rvec& s,
        jiba::rvec &Grad, const jiba::rvec &model, double &misfit, double *stp,
        int itnmax, double ftol, double xtol, double gtol, double stpmax,
        double stpmin);
    //! A backtracking line search taken from Opt++
    int backtrack(jiba::ObjectiveFunction* nlp, jiba::rvec& search_dir,
        jiba::rvec& grad, const jiba::rvec &model, double &misfit, double *stp,
        int itnmax, double ftol, double stpmax, double stpmin);

  }
#endif /* MCSRCH_H_ */
