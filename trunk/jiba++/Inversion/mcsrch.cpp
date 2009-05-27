//------------------------------------------------------------------------
// Copyright (C) 1993:
// J.C. Meza
// Sandia National Laboratories
// meza@california.sandia.gov
//
// modifications by Max Moorkamp 2009
//------------------------------------------------------------------------

#include "mcsrch.h"
#include <iostream>
namespace OPTPP
  {
    using std::min;
    using std::max;

    int mcstep(double *stx, double *fx, double *dx, double *sty, double *fy,
        double *dy, double *stp, double fp, double dp, bool *brackt,
        double stpmin, double stpmax, int *info);

    int mcsrch(jiba::ObjectiveFunction* nlp, const jiba::rvec& s,
        jiba::rvec &Grad,
        const jiba::rvec &model, double misfit, double *stp, int itnmax,
        double ftol, double xtol, double gtol, double stpmax, double stpmin)
      {

        /****************************************************************************
         *   subroutine mcsrch
         *   Purpose
         *   find a step which satisfies
         *   a sufficient decrease condition and a curvature condition.
         *
         *   at each stage the subroutine updates an interval of
         *   uncertainty with endpoints stx and sty. the interval of
         *   uncertainty is initially chosen so that it contains a
         *   minimizer of the modified function
         *        f(x+stp*s) - f(x) - ftol*stp*(gradf(x)'s).
         *
         *   if a step is obtained for which the modified function
         *   has a nonpositive function value and nonnegative derivative,
         *   then the interval of uncertainty is chosen so that it
         *   contains a minimizer of f(x+stp*s).
         *   the algorithm is designed to find a step which satisfies
         *   the sufficient decrease condition
         *         f(x+stp*s) .le. f(x) + ftol*stp*(gradf(x)'s),
         *   and the curvature condition
         *         abs(gradf(x+stp*s)'s)) .le. gtol*abs(gradf(x)'s).
         *   if ftol is less than gtol and if, for example, the function
         *   is bounded below, then there is always a step which satisfies
         *   both conditions. if no step can be found which satisfies both
         *   conditions, then the algorithm usually stops when rounding
         *   errors prevent further progress. in this case stp only
         *   satisfies the sufficient decrease condition.
         *
         *   Parameters
         *     s is an input array of length n which specifies the
         *       search direction.
         *     stp is a nonnegative variable. on input stp contains an
         *       initial estimate of a satisfactory step. on output
         *       stp contains the final estimate.
         *     ftol and gtol are nonnegative input variables.
         *       termination occurs when the sufficient decrease
         *       condition and the directional derivative condition are
         *       satisfied.
         *       ftol should be smaller than 5.e-1
         *       suggested value = 1.e-4 for newton methods
         *                       = 1.e-1 for more exact line searches
         *       Default Value = 1.e-4
         *
         *       gtol should be greater than 1.e-4
         *       Default Value = 0.9
         *     xtol is a nonnegative input variable. termination occurs
         *       when the relative width of the interval of uncertainty
         *       is at most xtol.
         *       Default Value = 2.2e-16
         *     stpmin and stpmax are nonnegative input variables which
         *       specify lower and upper bounds for the step. (in this reverse
         *       communication implementatin they are defined in a common
         *       statement).
         *       stpmin Default Value = 1.e-9
         *       stpmax Default Value = 1.e3
         *     info is an integer output variable set as follows:
         *       info =-1 improper input parameters.
         *       info = 1  the sufficient decrease condition and the
         *                 directional derivative condition hold.
         *       info =-2  relative width of the interval of uncertainty
         *                 is at most xtol.
         *       info =-4  the step is at the lower bound stpmin.
         *       info =-5  the step is at the upper bound stpmax.
         *       info =-6  rounding errors prevent further progress.
         *                 there may not be a step which satisfies the
         *                 sufficient decrease and curvature conditions.
         *                 tolerances may be too small.
         *   subprograms called
         *     mcstep
         *     fortran-supplied...abs,max,min
         *   argonne national laboratory. minpack project. june 1983
         *   jorge j. more', david j. thuente
         *
         *
         *   Recoded in C++ by Juan Meza December 1992
         *
         *
         *****************************************************************************/

        /* initialized data */

        static const double half = .5;
        static const double p66 = .66;
        static const double xtrapf = 4.;
        static const double zero = 0.;

        /* local variables */
        double dgxm, dgym;
        int j, info, infoc;
        double finit, width, stmin, stmax;
        bool stage1;
        double width1, ftest1, dg, fm, fx, fy;
        bool brackt;
        double dginit, dgtest;
        double dgm, dgx, dgy, fxm, fym, stx, sty;

        int siter;
        //int    maxiter = itnmax;
        int maxiter = 10;
        int n = s.size();

        double fvalue;
        jiba::rvec xc(n);

        infoc = 1;

        /*   check the input parameters for errors. */

        if (n <= 0 || *stp <= zero || ftol < zero || gtol < zero || xtol < zero
            || stpmin < zero || stpmax < stpmin)
          {
            infoc = -1;
            return infoc;
          }

        /* compute the initial gradient in the search direction */
        /* and check that s is a descent direction. */

        dginit = zero;

        for (j = 0; j < n; ++j)
          {
            dginit += Grad(j) * s(j);
          }

        if (dginit >= zero)
          {
            std::cout
                << "\nmcsrch: Initial search direction not a descent direction\n";
            return -1;
          }

        /* initialize local variables. */

        brackt = false;
        stage1 = true;
        finit = misfit;
        dgtest = ftol * dginit;
        width = stpmax - stpmin;
        width1 = width / half;

        /* the variables stx, fx, dgx contain the values of the step, */
        /* function, and directional derivative at the best step. */
        /* the variables sty, fy, dgy contain the value of the step, */
        /* function, and derivative at the other endpoint of */
        /* the interval of uncertainty. */
        /* the variables stp, f, dg contain the values of the step, */
        /* function, and derivative at the current step. */

        stx = zero;
        fx = finit;
        dgx = dginit;
        sty = zero;
        fy = finit;
        dgy = dginit;

        siter = 0;

        /* start of iteration. */

        while (siter < maxiter)
          {
            siter++;

            /*  set the minimum and maximum steps to correspond */
            /*  to the present interval of uncertainty. */

            if (brackt)
              {
                stmin = std::min(stx, sty);
                stmax = std::max(stx, sty);
              }
            else
              {
                stmin = stx;
                stmax = *stp + xtrapf * (*stp - stx);
              }

            /*    force the step to be within the bounds stpmax and stpmin. */

            *stp = std::max(*stp, stpmin);
            *stp = std::min(*stp, stpmax);

            /*    if an unusual termination is to occur then let */
            /*    stp be the lowest point obtained so far. */

            if (brackt && (*stp <= stmin || *stp >= stmax) || infoc == 0
                || brackt && stmax - stmin <= xtol * stmax)
              {
                *stp = stx;
              }

            /*    evaluate the function and gradient at stp */
            /*    and compute the directional derivative. */

            xc = model + s * (*stp);

            fvalue = nlp->CalcMisfit(xc);
            Grad = nlp->CalcGradient();

            info = 0;
            dg = zero;
            for (j = 0; j < n; ++j)
              {
                dg += Grad(j) * s(j);
              }
            ftest1 = finit + *stp * dgtest;

            /*    test for convergence. */

            if (brackt && (*stp <= stmin || *stp >= stmax) || infoc == 0)
              {
                info = -6; // CPJW 12/10/2003 original stmt info = 6
              }
            if (*stp == stpmax && fvalue <= ftest1 && dg <= dgtest)
              {
                info = -5; // CPJW 12/10/2003 original stmt info = 5
              }
            if (*stp == stpmin && (fvalue > ftest1 || dg >= dgtest))
              {
                info = -4; // CPJW 12/10/2003 original stmt info = 4
              }
            if (brackt && stmax - stmin <= xtol * stmax)
              {
                info = -2; // CPJW 12/10/2003 original stmt info = 2
              }
            if (fvalue <= ftest1 && fabs(dg) <= gtol * (-dginit))
              {
                info = 1;
              }

            /*    check for termination. */

            if (info < 0)
              {
                return (-1);
              }
            if (info != 0)
              {
                if (siter == 1)
                  return (0);
                else
                  return (1);
              }

            /*  in the first stage we seek a step for which the modified */
            /*  function has a nonpositive value and nonnegative derivative. */

            if (stage1 && fvalue <= ftest1 && dg >= min(ftol, gtol) * dginit)
              {
                stage1 = false;
              }

            /*  a modified function is used to predict the step only if */
            /*  we have not obtained a step for which the modified */
            /*  function has a nonpositive function value and nonnegative */
            /*  derivative, and if a lower function value has been */
            /*  obtained but the decrease is not sufficient. */

            if (stage1 && fvalue <= fx && fvalue > ftest1)
              {

                /* define the modified function and derivative values. */

                fm = fvalue - *stp * dgtest;
                fxm = fx - stx * dgtest;
                fym = fy - sty * dgtest;
                dgm = dg - dgtest;
                dgxm = dgx - dgtest;
                dgym = dgy - dgtest;

                /* call cstep to update the interval of uncertainty */
                /* and to compute the new step. */

                mcstep(&stx, &fxm, &dgxm, &sty, &fym, &dgym, stp, fm, dgm,
                    &brackt, stmin, stmax, &infoc);

                /* reset the function and gradient values for f. */

                fx = fxm + stx * dgtest;
                fy = fym + sty * dgtest;
                dgx = dgxm + dgtest;
                dgy = dgym + dgtest;
              }
            else
              {

                /* call mcstep to update the interval of uncertainty */
                /* and to compute the new step. */

                mcstep(&stx, &fx, &dgx, &sty, &fy, &dgy, stp, fvalue, dg,
                    &brackt, stmin, stmax, &infoc);
              }

            /*  force a sufficient decrease in the size of the */
            /*  interval of uncertainty. */

            if (brackt)
              {
                if (fabs(sty - stx) >= p66 * width1)
                  {
                    *stp = stx + (sty - stx) / 2.;
                  }
                width1 = width;
                width = fabs(sty - stx);
              }

          } /*  end of iteration. */

        return (-1); // too many iterations;

      } /* last line of subroutine mcsrch. */
    /****************************************************************************/
    int mcstep(double *stx, double *fx, double *dx, double *sty, double *fy,
        double *dy, double *stp, double fp, double dp, bool *brackt,
        double stpmin, double stpmax, int *info)
      {
        /* system generated locals */
        double d_1, d_2, d_3;

        /* local variables */
        static double sgnd, stpc, stpf, stpq, p, q, gamma, r, s, theta;
        static bool bound;

        /******************************************************************************
         *     subroutine mcstep

         *     the purpose of mcstep is to compute a safeguarded step for
         *     a linesearch and to update an interval of uncertainty for
         *     a minimizer of the function.

         *     the parameter stx contains the step with the least function
         *     value. the parameter stp contains the current step. it is
         *     assumed that the derivative at stx is negative in the
         *     direction of the step. if brackt is set true then a
         *     minimizer has been bracketed in an interval of uncertainty
         *     with endpoints stx and sty.

         *     the subroutine statement is

         *       subroutine mcstep(stx,fx,dx,sty,fy,dy,stp,fp,dp,brackt,
         *                        stpmin,stpmax,info)

         *     where

         *       stx, fx, and dx are variables which specify the step,
         *         the function, and the derivative at the best step obtained
         *         so far. the derivative must be negative in the direction
         *         of the step, that is, dx and stp-stx must have opposite
         *         signs. on output these parameters are updated appropriately.


         *       sty, fy, and dy are variables which specify the step,
         *         the function, and the derivative at the other endpoint of
         *         the interval of uncertainty. on output these parameters are
         *         updated appropriately.

         *       stp, fp, and dp are variables which specify the step,
         *         the function, and the derivative at the current step.
         *         if brackt is set true then on input stp must be
         *         between stx and sty. on output stp is set to the new step.

         *       brackt is a long int variable which specifies if a minimizer
         *         has been bracketed. if the minimizer has not been bracketed
         *         then on input brackt must be set false. if the minimizer
         *         is bracketed then on output brackt is set true.

         *       stpmin and stpmax are input variables which specify lower
         *         and upper bounds for the step.

         *       info is an int output variable set as follows:
         *         if info = 1,2,3,4,5, then the step has been computed
         *         according to one of the five cases below. otherwise
         *         info = 0, and this indicates improper input parameters.

         *     subprograms called

         *       fortran-supplied ... abs,max,min,sqrt

         *     argonne national laboratory. minpack project. june 1983
         *     jorge j. more', david j. thuente

         ******************************************************************************/

        *info = 0;

        /* check the input parameters for errors. */

        if (*brackt && (*stp <= min(*stx, *sty) || *stp >= max(*stx, *sty))
            || *dx * (*stp - *stx) >= 0. || stpmax < stpmin)
          {
            return 0;
          }

        /* determine if the derivatives have opposite sign. */

        sgnd = dp * (*dx / fabs(*dx));

        /* first case. a higher function value. */
        /* the minimum is bracketed. if the cubic step is closer */
        /* to stx than the quadratic step, the cubic step is taken, */
        /* else the average of the cubic and quadratic steps is taken. */

        if (fp > *fx)
          {
            *info = 1;
            bound = true;
            theta = (*fx - fp) * 3 / (*stp - *stx) + *dx + dp;
            /* computing max */
            d_1 = fabs(theta), d_2 = fabs(*dx), d_1 = max(d_2, d_1), d_2
                = fabs(dp);
            s = max(d_2, d_1);
            /* computing 2nd power */
            d_1 = theta / s;
            gamma = s * sqrt(d_1 * d_1 - *dx / s * (dp / s));
            if (*stp < *stx)
              {
                gamma = -gamma;
              }
            p = gamma - *dx + theta;
            q = gamma - *dx + gamma + dp;
            r = p / q;
            stpc = *stx + r * (*stp - *stx);
            stpq = *stx + *dx / ((*fx - fp) / (*stp - *stx) + *dx) / 2 * (*stp
                - *stx);
            if ((d_1 = stpc - *stx, fabs(d_1)) < (d_2 = stpq - *stx, fabs(d_2)))
              {
                stpf = stpc;
              }
            else
              {
                stpf = stpc + (stpq - stpc) / 2;
              }
            *brackt = true;

            /* second case. a lower function value and derivatives of */
            /* opposite sign. the minimum is bracketed. if the cubic */
            /* step is closer to stx than the quadratic (secant) step, */
            /* the cubic step is taken, else the quadratic step is taken. */

          }
        else if (sgnd < 0.)
          {
            *info = 2;
            bound = false;
            theta = (*fx - fp) * 3 / (*stp - *stx) + *dx + dp;
            /* computing max */
            d_1 = fabs(theta), d_2 = fabs(*dx), d_1 = max(d_2, d_1), d_2
                = fabs(dp);
            s = max(d_2, d_1);
            /* computing 2nd power */
            d_1 = theta / s;
            gamma = s * sqrt(d_1 * d_1 - *dx / s * (dp / s));
            if (*stp > *stx)
              {
                gamma = -gamma;
              }
            p = gamma - dp + theta;
            q = gamma - dp + gamma + *dx;
            r = p / q;
            stpc = *stp + r * (*stx - *stp);
            stpq = *stp + dp / (dp - *dx) * (*stx - *stp);
            if ((d_1 = stpc - *stp, fabs(d_1)) > (d_2 = stpq - *stp, fabs(d_2)))
              {
                stpf = stpc;
              }
            else
              {
                stpf = stpq;
              }
            *brackt = true;

            /* third case. a lower function value, derivatives of the */
            /* same sign, and the magnitude of the derivative decreases. */
            /* the cubic step is only used if the cubic tends to infinity */
            /* in the direction of the step or if the minimum of the cubic */
            /* is beyond stp. otherwise the cubic step is defined to be */
            /* either stpmin or stpmax. the quadratic (secant) step is also */

            /* computed and if the minimum is bracketed then the the step */
            /* closest to stx is taken, else the step farthest away is taken.
             */

          }
        else if (fabs(dp) < fabs(*dx))
          {
            *info = 3;
            bound = true;
            theta = (*fx - fp) * 3 / (*stp - *stx) + *dx + dp;
            /* computing max */
            d_1 = fabs(theta), d_2 = fabs(*dx), d_1 = max(d_2, d_1), d_2
                = fabs(dp);
            s = max(d_2, d_1);

            /* the case gamma = 0 only arises if the cubic does not tend */
            /* to infinity in the direction of the step. */

            /* computing max */
            /* computing 2nd power */
            d_3 = theta / s;
            d_1 = 0., d_2 = d_3 * d_3 - *dx / s * (dp / s);
            gamma = s * sqrt((max(d_2, d_1)));
            if (*stp > *stx)
              {
                gamma = -gamma;
              }
            p = gamma - dp + theta;
            q = gamma + (*dx - dp) + gamma;
            r = p / q;
            if (r < 0. && gamma != 0.)
              {
                stpc = *stp + r * (*stx - *stp);
              }
            else if (*stp > *stx)
              {
                stpc = stpmax;
              }
            else
              {
                stpc = stpmin;
              }
            stpq = *stp + dp / (dp - *dx) * (*stx - *stp);
            if (*brackt)
              {
                if ((d_1 = *stp - stpc, fabs(d_1)) < (d_2 = *stp - stpq, fabs(
                    d_2)))
                  {
                    stpf = stpc;
                  }
                else
                  {
                    stpf = stpq;
                  }
              }
            else
              {
                if ((d_1 = *stp - stpc, fabs(d_1)) > (d_2 = *stp - stpq, fabs(
                    d_2)))
                  {
                    stpf = stpc;
                  }
                else
                  {
                    stpf = stpq;
                  }
              }

            /* fourth case. a lower function value, derivatives of the */
            /* same sign, and the magnitude of the derivative does */
            /* not decrease. if the minimum is not bracketed, the step */
            /* is either stpmin or stpmax, else the cubic step is taken. */

          }
        else
          {
            *info = 4;
            bound = false;
            if (*brackt)
              {
                theta = (fp - *fy) * 3 / (*sty - *stp) + *dy + dp;
                /* computing max */
                d_1 = fabs(theta), d_2 = fabs(*dy), d_1 = max(d_2, d_1), d_2
                    = fabs(dp);
                s = max(d_2, d_1);
                /* computing 2nd power */
                d_1 = theta / s;
                gamma = s * sqrt(d_1 * d_1 - *dy / s * (dp / s));
                if (*stp > *sty)
                  {
                    gamma = -gamma;
                  }
                p = gamma - dp + theta;
                q = gamma - dp + gamma + *dy;
                r = p / q;
                stpc = *stp + r * (*sty - *stp);
                stpf = stpc;
              }
            else if (*stp > *stx)
              {
                stpf = stpmax;
              }
            else
              {
                stpf = stpmin;
              }
          }

        /* update the interval of uncertainty. this update does not */
        /* depend on the new step or the case analysis above. */

        if (fp > *fx)
          {
            *sty = *stp;
            *fy = fp;
            *dy = dp;
          }
        else
          {
            if (sgnd < 0.)
              {
                *sty = *stx;
                *fy = *fx;
                *dy = *dx;
              }
            *stx = *stp;
            *fx = fp;
            *dx = dp;
          }

        /*  compute the new step and safeguard it. */

        stpf = min(stpmax, stpf);
        stpf = max(stpmin, stpf);
        *stp = stpf;
        if (*brackt && bound)
          {
            if (*sty > *stx)
              {
                /* computing max */
                d_1 = *stx + (*sty - *stx) * (float) .66;
                *stp = min(*stp, d_1);
              }
            else
              {
                /* computing max */
                d_1 = *stx + (*sty - *stx) * (float) .66;
                *stp = max(*stp, d_1);
              }
          }
        return 0;

      }

  } // namespace OPTPP
