/*
 * Entropy.cpp
 *
 *  Created on: 20 Oct 2020
 *      Author: moorkamp
 */

#include "Entropy.h"

namespace jif3D
  {
    //Calculate the Shannon Entropy for a vector of values x
    //Typically these x come from a histogram or similar
    // def shan_entropy(c):
    //     c_normalized = c / float(np.sum(c))
    //     c_normalized = c_normalized[np.nonzero(c_normalized)]
    //     H = -sum(c_normalized* np.log(c_normalized))
    //     return H
    double shan_entropy(const jif3D::rvec &x)
      {
        double H = 0;
        const size_t nx = x.size();
        const double sum = ublas::sum(x);
#pragma omp parallel for reduction(-:H)
        for (size_t i = 0; i < nx; ++i)
          {
            // lim x->0 of x log(x) = 0 however for small values of x
            // log(x) overflows, we avoid this by arbitrarily cutting of
            // at 1e-100, experiments show that this avoids all problems
            if (x(i) > 1e-100)
              {
                const double val = x(i) / sum;
                H -= val * std::log(val);
              }
          }
        return H;
      }

    //calculate the derivative of the shannon entropy
    jif3D::rvec diff_shan_entropy(const jif3D::rvec &x)
      {
        //related output from maxima
        // -->  diff(x*log(x),x,1);
        //(%o1) log(x)+1
        // --> t:x/(a+b+x);
        // (t) x/(x+b+a)
        //  -->  diff(t*log(t),x,1);
        // (%o14)  log(x/(x+b+a))/(x+b+a)-(x*log(x/(x+b+a)))/(x+b+a)^2+1/(x+b+a)-x/(x+b+a)^2

        const size_t nx = x.size();
        const double sum = ublas::sum(x);
        jif3D::rvec result(nx, 0.0);
        for (size_t i = 0; i < nx; ++i)
          {
            //see explanation in shannon_entropy for this cut off
            if (x(i) > 1e-100)
              {
                const double val = x(i) / sum;
                //we have to do a chain rule here and consider
                //that for the derivative from the sum we get terms for all x_i
                //as the the variable sum depends on x
                result(i) -= (std::log(val) + 1.0) / sum;
                for (size_t j = 0; j < nx; ++j)
                  {
                    result(j) += (x(i) * std::log(val) + x(i)) / (sum * sum);
                  }
              }
          }
        return result;
      }
  }
