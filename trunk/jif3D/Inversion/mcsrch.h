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
    int mcsrch(jif3D::ObjectiveFunction* nlp, const jif3D::rvec& s,
        jif3D::rvec &Grad, const jif3D::rvec &model, double &misfit, double *stp,
        int itnmax, double ftol, double xtol, double gtol, double stpmax,
        double stpmin, bool Verbose = false);
    //! A backtracking line search taken from Opt++
    int backtrack(jif3D::ObjectiveFunction* nlp, jif3D::rvec& search_dir,
        jif3D::rvec& grad, const jif3D::rvec &model, double &misfit, double *stp,
        int itnmax, double ftol, double stpmax, double stpmin);

  }
#endif /* MCSRCH_H_ */
