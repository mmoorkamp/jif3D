/*
 * fittingtest.cpp
 *
 *  Created on: Aug 26, 2008
 *      Author: mmoorkamp
 */

#include <iostream>
#include "../Inversion/LinearInversion.h"
#include "../Global/Interpolate.h"

double func(const double z, const double z0)
  {
    return pow(z - z0, -2);
  }

double deriv(const double z, const double z0)
  {
    return -2.0 * pow(z - z0,-3);
  }

double misfit(const jiba::rvec &observed, const jiba::rvec &zvalues,
    const jiba::rvec &Model)
  {
    double fit = 0;
    for (size_t i = 0; i < observed.size(); ++i)
      {
        fit += pow(observed(i) - func(zvalues(i), Model(0)), 2);
      }
    return fit;
  }

int main()
  {
    const size_t ndata = 20;
    const double evalthresh = 1e-3;
    jiba::rvec observed(ndata), zvalues(ndata), error(ndata), delta(ndata);
    const double z0 = 350;
    const double zstart = 250;
    for (size_t i = 0; i < ndata; ++i)
      {
        zvalues(i) = i * 10;
        observed(i) = func(zvalues(i),z0);
        error(i) = 1.0;
      }
    jiba::rmat sens(ndata, 1);
    jiba::rvec weights(1);
    jiba::rvec InvModel(1), DeltaModel(1);
    weights(0) = 1.0;
    InvModel(0) = zstart;
    const size_t iterations = 10;

    for (size_t i = 0; i < iterations; ++i)
      {
        for (size_t j = 0; j < ndata; ++j)
          {
            sens(j, 0) = deriv(zvalues(j), InvModel(0));
            delta(j) = observed(j) - func(zvalues(j),InvModel(0));
          }
        jiba::ModelSpaceInversion()(sens, delta, weights, error,1e-5, 1e-10,
            DeltaModel);
        double stepsize = 1.0;
        //jiba::rvec TryModel1 = InvModel - stepsize * DeltaModel;
       // double m1 = misfit(observed, zvalues, TryModel1);
        //jiba::rvec TryModel2 = InvModel - 0.1 * stepsize * DeltaModel;
        //double m2 = misfit(observed, zvalues, TryModel2);
        //jiba::rvec TryModel3 = InvModel -  10.0 * stepsize * DeltaModel;
        //double m3 = misfit(observed, zvalues, TryModel3);
        //stepsize = jiba::QuadraticInterpolation(stepsize,m1,stepsize*0.1,m2,stepsize*10,m3);
        InvModel -= stepsize * DeltaModel;
        std::cout << InvModel << " " << DeltaModel << " " << stepsize
            << std::endl;
      }
  }
