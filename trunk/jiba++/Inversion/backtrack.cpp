//------------------------------------------------------------------------
// Copyright (C) 1993:
// J.C. Meza
// Sandia National Laboratories
// meza@california.sandia.gov
//------------------------------------------------------------------------


#ifdef HAVE_STD
#include <cmath>
#else
#include <math.h>
#endif
#include "ObjectiveFunction.h"
#include "mcsrch.h"
//------------------------------------------------------------------------
// external subroutines referenced by this module
//------------------------------------------------------------------------

#if !(defined(__GNUC__) && __GNUC__ >= 3)
extern "C"
  {
    double copysign(double, double);
  }
#else
#ifdef CYGWIN
extern "C"
  {
    extern double copysign _PARAMS((double, double));
  }
#endif
#endif

namespace OPTPP
  {

    int backtrack(jiba::ObjectiveFunction* nlp, jiba::rvec& search_dir,
        jiba::rvec& grad, const jiba::rvec &model, double &misfit, double *stp,
        int itnmax, double ftol, double stpmax, double stpmin)
      {
        /* Local variables */
        double disc;
        double a, b;
        int i;
        double fplus, t1, t2, t3, tlmbda, minlambda;
        double scl, rellength, sln, initslope;

        double lambda = 1.;
        double plmbda = -1.0;
        double pfplus = -1.0;
        double one;

        int n = search_dir.size();
        jiba::rvec xplus(n);
        double fx, tmp1, tmp2;
        int iter;

        sln = ublas::norm_2(search_dir);

        if (sln >= stpmax && sln != 0.0)
          { // STEP LONGER THAN MAXIMUM ALLOWED
            scl = stpmax / sln;
            search_dir *= scl;
            sln = stpmax;
          }

        initslope = ublas::inner_prod(grad, search_dir);
        if (initslope >= 0.0)
          {
            search_dir = -grad;
            initslope = -ublas::inner_prod(grad, grad);
            std::cout << "Not a descent direction !" << std::endl;
          }

        rellength = 0.;

        for (i = 0; i < n; ++i)
          {
            tmp1 = fabs(search_dir(i));
            tmp2 = fabs(model(i));
            rellength = std::max(rellength, tmp1) / tmp2;
          }

        minlambda = stpmin; // / rellength;

        /* CHECK IF NEW ITERATE SATISFACTORY.  GENERATE NEW LAMBDA IF NECESSARY. */
        iter = 0;
        //const  int MAXITER = itnmax;
        // std::cout << "Search Dir: " << search_dir << std::endl;
        while (iter < itnmax)
          {

            iter++;
            xplus = model + search_dir * (lambda);

            fplus = nlp->CalcMisfit(xplus);
            //std::cout << "Lambda: " << lambda << " Misfit: " << fplus << " Model: " << xplus <<  std::endl;
            // Is this an acceptable step ?

            if (fplus <= misfit + initslope * ftol * lambda)
              {
                // Yes !
                *stp = lambda;
                misfit = fplus;
                grad = nlp->CalcGradient(xplus);

                if (iter == 1)
                  return (0); // Newton Step
                else
                  return (1); // Backtrack Step
              }

            // Insufficient decrease. Try to backtrack

            if (lambda < minlambda)
              { // Step size smaller than min allowed
                *stp = lambda;
                misfit = nlp->CalcMisfit(model);
                grad = nlp->CalcGradient(model);
                std::cout << "Step size smaller than min allowed" << std::endl;
                return (-1); // Error
              }
            else
              { // Compute new lambda

                if (iter == 1)
                  {/* FIRST BACKTRACK: QUADRATIC FIT */

                    tlmbda = -initslope / ((fplus - fx - initslope) * 2.);
                  }
                else
                  {/* ALL SUBSEQUENT BACKTRACKS: CUBIC FIT */

                    t1 = fplus - fx - lambda * initslope;
                    t2 = pfplus - fx - plmbda * initslope;
                    t3 = 1. / (lambda - plmbda);
                    a = t3 * (t1 / (lambda * lambda) - t2 / (plmbda * plmbda));
                    b = t3 * (t2 * lambda / (plmbda * plmbda) - t1 * plmbda
                        / (lambda * lambda));
                    disc = b * b - a * 3. * initslope;
                    if (disc > b * b)
                      {
                        // ONLY ONE POSITIVE CRITICAL POINT, MUST BE MINIMUM
                        one = copysign(1.0, a);
                        tlmbda = (-b + one * sqrt(disc)) / (a * 3.);
                      }
                    else
                      {
                        // BOTH CRITICAL POINTS POSITIVE, FIRST IS MINIMUM
                        one = copysign(1.0, a);
                        tlmbda = (-(float) b - one * sqrt(disc)) / (a * 3.);
                      }
                    if (tlmbda > lambda * .5)
                      tlmbda = lambda * .5;
                  }

                plmbda = lambda;
                pfplus = fplus;
                if (tlmbda < .1 * lambda)
                  lambda *= .1;
                else
                  lambda = tlmbda;
              }
          }
        std::cout << "Too many iterations" << std::endl;
        misfit = nlp ->CalcMisfit(model);
        grad = nlp->CalcGradient(model);
        return (-1); // Too many iterations
      }

  } // namespace OPTPP
