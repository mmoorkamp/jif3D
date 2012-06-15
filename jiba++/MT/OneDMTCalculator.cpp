//============================================================================
// Name        : OneDMTCalculator.cpp
// Author      : 13 Jun 2012
// Version     : 
// Copyright   : 2012, mm489
//============================================================================

#include <boost/math/special_functions/fpclassify.hpp>
#include "OneDMTCalculator.h"
#include "MTEquations.h"

namespace jiba
  {

    OneDMTCalculator::OneDMTCalculator()
      {
        // TODO Auto-generated constructor stub

      }

    OneDMTCalculator::~OneDMTCalculator()
      {
        // TODO Auto-generated destructor stub
      }

    jiba::rvec OneDMTCalculator::Calculate(const ModelType &Model)
      {
        const size_t nfreq = Model.GetFrequencies().size();
        const size_t nlayers = Model.GetBackgroundThicknesses().size();
        const double PI = std::acos(-1.0);
        const std::complex<double> I(0.0, 1.0);
        std::complex<double> omegamu;
        std::complex<double> kcurr, klow;
        std::complex<double> xi;
        std::vector<std::complex<double> > alpha(nlayers, 0);
        std::complex<double> adm;
        double d;
        double sigmacurr, sigmalow;

        jiba::rvec result(nfreq * 2);
        Z.resize(nfreq);
        Gamma.resize(nlayers, nfreq);
        gammakj.resize(nlayers, nfreq);
        gammaj.resize(nlayers, nfreq);
        for (size_t i = 0; i < nfreq; ++i)
          {
            const double omega = 2. * PI * Model.GetFrequencies().at(i);
            omegamu = -I * 4e-7 * PI * omega;
            fill(alpha.begin(), alpha.end(), 0);
            sigmacurr = Model.GetBackgroundConductivities().back();
            kcurr = sqrt(omegamu * sigmacurr);
            gammakj(nlayers - 1, i) = 0.5 / kcurr;
            for (int layerindex = nlayers - 2; layerindex >= 0; --layerindex)
              {
                sigmalow = Model.GetBackgroundConductivities().at(layerindex + 1);
                sigmacurr = Model.GetBackgroundConductivities().at(layerindex);
                kcurr = sqrt(omegamu * sigmacurr);
                klow = sqrt(omegamu * sigmalow);
                if (kcurr.real() < 0.0)
                  {
                    kcurr *= -1.0;
                    klow *= -1.0;
                  }
                xi = omegamu * (sigmacurr - sigmalow) / std::pow(kcurr + klow, 2);
                alpha.at(layerindex + 1) = (xi + alpha.at(layerindex + 1))
                    / (1. + xi * alpha.at(layerindex + 1));
                d = Model.GetBackgroundThicknesses().at(layerindex);
                std::complex<double> dk = kcurr * d;
                std::complex<double> xi1 = exp(-dk);
                alpha.at(layerindex) = alpha.at(layerindex + 1) * xi1 * xi1;
                gammaj(layerindex, i) = (1. + alpha.at(layerindex + 1))
                    / (1. / xi1 + xi1 * alpha.at(layerindex + 1));
                std::complex<double> mu = alpha.at(layerindex);
                std::complex<double> c = 1. + mu;
                gammakj(layerindex, i) = d
                    * (2.0 * mu / (c * c)
                        + 0.5
                            * (std::pow(mu / (c * xi1), 2) - std::pow(xi1 / c, 2)
                                + (1. - mu) / c) / dk);
              }
            adm = kcurr / (-I * omega) * ((1. - alpha.at(0)) / (1. + alpha.at(0)));
            Z(i) = jiba::mag_mu * conj(1. / adm);
            result(2 * i) = Z(i).real();
            result(2 * i + 1) = Z(i).imag();
          }
        return result;
      }

    jiba::rvec OneDMTCalculator::LQDerivative(const ModelType &Model, const rvec &Misfit)
      {
        const size_t nfreq = Model.GetFrequencies().size();
        const size_t nlayers = Model.GetBackgroundThicknesses().size();
        jiba::rvec result(nlayers, 0.0);
        for (size_t i = 0; i < nfreq; ++i)
          {
            Gamma(0, i) = gammakj(0, i);
            for (size_t layerindex = 0; layerindex < nlayers - 1; ++layerindex)
              {
                Gamma(layerindex + 1, i) = std::pow(gammaj(layerindex, i), 2)
                    * gammakj(layerindex + 1, i) / gammakj(layerindex, i)
                    * Gamma(layerindex, i);
              }
            for (size_t j = 0; j < nlayers; ++j)
              {
                result(j) += std::real(
                    Gamma(j, i) * Z(i) * Z(i)
                        * std::complex<double>(Misfit(2 * i), Misfit(2 * i + 1)));
                if (boost::math::isnan(result(j)))
                  {
                    std::cout << j << " " << i << " " << Model.GetFrequencies().at(i) << " " << Gamma(j, i) << " " << Z(i)
                        << std::endl;
                  }
              }
          }
        return 2.0 * result;
      }

  } /* namespace jiba */
