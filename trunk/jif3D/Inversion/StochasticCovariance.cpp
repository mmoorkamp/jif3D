/*
 * StochasticCovariance.cpp
 *
 *  Created on: 15 Jun 2018
 *      Author: mm489
 */

#include "StochasticCovariance.h"
#include "../Global/NumUtil.h"
#include "../Global/kissfft.hh"
#include <boost/math/special_functions/gamma.hpp>
#include <boost/math/special_functions/bessel.hpp>
#include <fstream>
namespace jif3D
  {

    StochasticCovariance::StochasticCovariance(size_t x, size_t y, size_t z, double ma,
        double mnu, double msigma) :
        nx(x), ny(y), nz(z), a(ma), nu(mnu), sigma(msigma), HaveInv(false)
      {
        double factor = sigma * sigma / (std::pow(2, nu - 1) * boost::math::tgamma(nu));

        const size_t nmod = nx * ny * nz;
        Cm.resize(nmod, nmod);
        Cmf.resize(nmod * 2);
        Cmv.resize(nmod * 2);
#pragma omp parallel for default(shared)
        for (size_t i = 0; i < nmod; ++i)
          {
            int xi, yi, zi;
            OffsetToIndex(i, xi, yi, zi);
            //double currx = Model.GetXCoordinates()[xi] + Model.GetXCellSizes()[xi] / 2.0;
            //double curry = Model.GetYCoordinates()[yi] + Model.GetYCellSizes()[yi] / 2.0;
            //double currz = Model.GetZCoordinates()[zi] + Model.GetZCellSizes()[zi] / 2.0;
            for (size_t j = 0; j < nmod; ++j)
              {
                int xj, yj, zj;
                OffsetToIndex(j, xj, yj, zj);
                //double x = Model.GetXCoordinates()[xj] + Model.GetXCellSizes()[xj] / 2.0;
                //double y = Model.GetYCoordinates()[yj] + Model.GetYCellSizes()[yj] / 2.0;
                //double z = Model.GetZCoordinates()[zj] + Model.GetZCellSizes()[zj] / 2.0;
                //double r = std::sqrt(
                //    jif3D::pow2(currx - x) + jif3D::pow2(curry - y)
                //        + jif3D::pow2(currz - z));
                double r = std::sqrt(
                    jif3D::pow2(xi - xj) + jif3D::pow2(yi - yj) + jif3D::pow2(zi - zj));
                const double ra = std::abs(r / a);
                Cm(i, j) =
                    (ra == 0.0) ?
                        1.0 :
                        factor * std::pow(ra, nu) * boost::math::cyl_bessel_k(nu, ra);

              }

          }
        /*Cmv(0) = Cm(0, 0);
        for (size_t i = 1; i < nmod; ++i)
          {
            Cmv(i) = Cm(i, 0);
            Cmv(2*nmod- i) = Cm(i, 0);
          }
        Cmv(nmod) = Cmv(nmod+1);
        kissfft<double> fft(nmod * 2, false);
        fft.transform(&Cmv.data()[0], &Cmf.data()[0]);*/
        // TODO Auto-generated constructor stub
        /*std::ofstream ovfile("cmv.out");
        for (auto val : Cmv)
          {
            ovfile << std::abs(val) << "\n";
          }

        std::ofstream offile("cmf.out");
        for (auto val : Cmf)
          {
            offile << std::abs(val) << "\n";
          }*/
      }

    jif3D::rvec StochasticCovariance::ApplyCovar(const jif3D::rvec &vector)
      {

        double factor = sigma * sigma / (std::pow(2, nu - 1) * boost::math::tgamma(nu));

        const size_t nmod = nx * ny * nz;
        jif3D::rvec result(vector.size(),0.0);
#pragma omp parallel for default(shared)
        for (size_t i = 0; i < nmod; ++i)
          {
            int xi, yi, zi;
            OffsetToIndex(i, xi, yi, zi);
            //double currx = Model.GetXCoordinates()[xi] + Model.GetXCellSizes()[xi] / 2.0;
            //double curry = Model.GetYCoordinates()[yi] + Model.GetYCellSizes()[yi] / 2.0;
            //double currz = Model.GetZCoordinates()[zi] + Model.GetZCellSizes()[zi] / 2.0;
            for (size_t j = 0; j < nmod; ++j)
              {
                int xj, yj, zj;
                OffsetToIndex(j, xj, yj, zj);
                //double x = Model.GetXCoordinates()[xj] + Model.GetXCellSizes()[xj] / 2.0;
                //double y = Model.GetYCoordinates()[yj] + Model.GetYCellSizes()[yj] / 2.0;
                //double z = Model.GetZCoordinates()[zj] + Model.GetZCellSizes()[zj] / 2.0;
                //double r = std::sqrt(
                //    jif3D::pow2(currx - x) + jif3D::pow2(curry - y)
                //        + jif3D::pow2(currz - z));
                double r = std::sqrt(
                    jif3D::pow2(xi - xj) + jif3D::pow2(yi - yj) + jif3D::pow2(zi - zj));
                const double ra = std::abs(r / a);
                double matelem =
                    (ra == 0.0) ?
                        1.0 :
                        factor * std::pow(ra, nu) * boost::math::cyl_bessel_k(nu, ra);
                 result(i) += matelem * vector(j);
              }

          }


/*        typedef Eigen::Matrix<double, Eigen::Dynamic, 1> VectorXi;
        typedef Eigen::Map<const VectorXi> MapTypeConst;
        typedef Eigen::Map<VectorXi> MapType;

        MapTypeConst v(&(vector.data()[0]), vector.size(), 1);
        MapType r(&(result.data()[0]), result.size(), 1);
        r = Cm * v;*/

        return result;
        /*const size_t nmod = vector.size() * 2;
        jif3D::cvec inv(nmod, 0.0), invf(nmod, 0.0);
        jif3D::rvec result(vector.size(), 0.0);
        std::copy(vector.begin(), vector.end(), inv.begin());
        kissfft<double> fft(nmod, false);
        fft.transform(&inv.data()[0], &invf.data()[0]);
        invf = ublas::element_prod(invf, Cmf);
        kissfft<double> ifft(nmod, true);
        ifft.transform(&invf.data()[0], &inv.data()[0]);
        for (size_t i = 0; i < vector.size(); ++i)
          {
            result(i) = 1.0 / double(nmod) * inv(i).real();
          }
        return result;*/
      }

    jif3D::rvec StochasticCovariance::ApplyInvCovar(const jif3D::rvec &vector)
      {
        // if (!HaveInv)
        //   {
        //     solver.compute(Cm);
        //  }
        //HaveInv = true;
        //const size_t nmod = vector.size();
        //jif3D::rvec result(nmod);
        //typedef Eigen::Matrix<double, Eigen::Dynamic, 1> VectorXi;
        //typedef Eigen::Map<const VectorXi> MapTypeConst;
        //typedef Eigen::Map<VectorXi> MapType;

        //MapTypeConst v(&(vector.data()[0]), nmod, 1);
        //MapType r(&(result.data()[0]), nmod, 1);
        //r = solver.solve(v);
        //return result;

        /*const size_t nmod = vector.size() * 2;
        jif3D::cvec inv(nmod, 0.0), invf(nmod, 0.0);
        jif3D::rvec result(vector.size(), 0.0);
        std::copy(vector.begin(), vector.end(), inv.begin());
        kissfft<double> fft(nmod, false);
        fft.transform(&inv.data()[0], &invf.data()[0]);
       double maxCmf = std::abs(*std::max_element(Cmf.begin(), Cmf.end(),
            [](std::complex<double> a, std::complex<double> b)
              { return std::abs(a) > std::abs(b);}));
        double threshold = 1e-6 * maxCmf;
        std::ofstream invfile("inv.out");
        for (size_t i = 0; i < nmod; ++i)
          {
            invfile << i << " " << std::abs(invf(i));
            if (std::abs(Cmf(i)) < threshold)
              {
                invf(i) = invf(i) / threshold;
              }
            else
              {
                invf(i) = invf(i) / Cmf(i);
              }
            invfile << " " << std::abs(invf(i)) << " " << std::abs(Cmf(i)) << "\n";
          }*/
        //invf = ublas::element_div(invf, Cmf);
        /*kissfft<double> ifft(nmod, true);
        std::fill_n(inv.begin(),nmod,0.0);
        ifft.transform(&invf.data()[0], &inv.data()[0]);
        for (size_t i = 0; i < vector.size(); ++i)
          {
            result(i) = 1.0 / double(nmod) * inv(i).real();
          }
        return result;*/
        return vector;
      }

    StochasticCovariance::~StochasticCovariance()
      {
        // TODO Auto-generated destructor stub
      }

  } /* namespace jif3D */
