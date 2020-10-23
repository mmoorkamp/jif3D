/*
 * EntropyRegularization.cpp
 *
 *  Created on: 20 Oct 2020
 *      Author: moorkamp
 */

#include "EntropyRegularization.h"
#include "../MI/Entropy.h"
#include <boost/algorithm/minmax_element.hpp>
#include <boost/math/distributions/normal.hpp>

namespace jif3D
  {
    void GaussHist(const jif3D::rvec &x, double xmin, double xmax, size_t nbins,
        jif3D::rvec &CountsX)
      {

        boost::math::normal norm;
        const size_t nparm = x.size();

        double xw = (xmax - xmin) / nbins;
        const double binwidth = xw / 2.0;

        CountsX.resize(nbins);

        std::fill(CountsX.begin(), CountsX.end(), 0.0);

        for (size_t i = 0; i < nparm; ++i)
          {
            //int xc = std::floor(std::min((x(i) - xmin) / xw, nbins - 1));
            //int yc = std::floor(std::min((y(i) - ymin) / yw, nbins - 1));
#pragma omp parallel for
            for (size_t j = 0; j < nbins; ++j)
              {
                const double distx = x(i) - (xmin + (j + 0.5) * xw);
                const double val = boost::math::pdf(norm,
                    std::sqrt(distx * distx) / binwidth);

                CountsX(j) += val;
              }
          }
      }

    jif3D::rvec diff_Gauss(const jif3D::rvec &x, double xmin, double xmax, size_t nbins,
        const jif3D::rvec &CountsX)
      {
        boost::math::normal norm;
        const size_t nparm = x.size();
        double xw = (xmax - xmin) / nbins;
        const double binwidth = xw / 2.0;

        jif3D::rvec Result(nparm, 0.0);
        jif3D::rvec dsx = diff_shan_entropy(CountsX);

#pragma omp parallel for
        for (size_t i = 0; i < nparm; ++i)
          {
            //int xc = std::floor(std::min((x(i) - xmin) / xw, nbins - 1));
            //int yc = std::floor(std::min((y(i) - ymin) / yw, nbins - 1));
            for (size_t j = 0; j < nbins; ++j)
              {
                const double distx = x(i) - (xmin + (j + 0.5) * xw);

                const double val = boost::math::pdf(norm, distx / binwidth) / (binwidth * binwidth);
                Result(i) -= distx * val * dsx(j);

              }
          }

        return Result;
      }

    void EntropyRegularization::ImplDataDifference(const jif3D::rvec &Model,
        jif3D::rvec &Diff)
      {
        const size_t nbins = 100;

        CountsX.resize(nbins);
        GaussHist(Model, xmin, xmax, nbins, CountsX);
        double H_X = shan_entropy(CountsX);

        Diff.resize(1);
        Diff(0) = std::sqrt(H_X);
      }

//! The gradient of the regularization with respect to the model parameters
    jif3D::rvec EntropyRegularization::ImplGradient(const jif3D::rvec &Model,
        const jif3D::rvec &Diff)
      {

        jif3D::rvec result = diff_Gauss(Model, xmin, xmax,  nbins, CountsX);
        return result;
      }
  }
