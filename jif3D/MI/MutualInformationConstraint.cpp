/*
 * MutualInformationConstraint.cpp
 *
 *  Created on: 2 Jul 2019
 *      Author: moorkamp
 */

#include "MutualInformationConstraint.h"
#include <omp.h>
#if BOOST_VERSION >= 107000
#include <boost/histogram.hpp>
#endif

#include <boost/algorithm/minmax_element.hpp>
#include <boost/math/distributions/normal.hpp>
#include <boost/math/special_functions/fpclassify.hpp>
#include <algorithm>
#include <cmath>
#include <fstream>
#include "../Global/NetCDFTools.h"
#include "../Global/NumUtil.h"

using netCDF::NcDim;
using netCDF::NcFile;
using netCDF::NcVar;
using netCDF::NcVarAtt;
namespace jif3D
  {

    void plothist(const std::string &filename, double xmin, double xmax, double ymin,
        double ymax, size_t nbins, const jif3D::rvec &CountXY)
      {
        NcFile DataFile(filename, NcFile::replace);
        double xw = (xmax - xmin) / nbins;
        double yw = (ymax - ymin) / nbins;
        std::vector<double> x(nbins), y(nbins);
        for (size_t i = 0; i < nbins; ++i)
          {
            x[i] = (i + 1) * xw;
            y[i] = (i + 1) * yw;
          }

        NcDim XDim = DataFile.addDim("x", nbins);
        NcVar XVar = DataFile.addVar("x", netCDF::ncDouble, XDim);
        XVar.putVar(std::vector<std::size_t>(
          { 0 }), std::vector<std::size_t>(
          { nbins }), x.data());

        NcDim YDim = DataFile.addDim("y", nbins);
        NcVar YVar = DataFile.addVar("y", netCDF::ncDouble, YDim);
        YVar.putVar(std::vector<std::size_t>(
          { 0 }), std::vector<std::size_t>(
          { nbins }), y.data());

        std::vector<NcDim> dims;
        dims.push_back(XDim);
        dims.push_back(YDim);

        NcVar ZVar = DataFile.addVar("z", netCDF::ncDouble, dims);
        ZVar.putVar(&CountXY(0));

      }
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

    double MIHist(const jif3D::rvec &x, const jif3D::rvec &y, double xmin, double xmax,
        double ymin, double ymax, size_t nbins, jif3D::rvec &CountsX,
        jif3D::rvec &CountsY, jif3D::rvec &CountsXY)
      {
#if BOOST_VERSION >= 107000
        using namespace boost::histogram;
        using namespace literals;
        auto hxy = make_histogram(axis::regular<>(nbins, xmin, xmax, "x"),
            axis::regular<>(nbins, ymin, ymax, "y"));
        const size_t nparm = x.size();
        assert(y.size() == nparm);
        for (size_t i = 0; i < nparm; ++i)
          {
            hxy(x(i), y(i));
          }
        auto hx = algorithm::project(hxy, 0_c);
        auto hy = algorithm::project(hxy, 1_c);

        CountsX.resize(nbins);
        CountsY.resize(nbins);
        CountsXY.resize(nbins * nbins);
        std::fill(CountsX.begin(), CountsX.end(), 0.0);
        std::fill(CountsY.begin(), CountsY.end(), 0.0);
        std::fill(CountsXY.begin(), CountsXY.end(), 0.0);
        for (auto&& xb : indexed(hx))
          {
            CountsX(xb.index()) = *xb;
          }
        for (auto&& yb : indexed(hy))
          {
            CountsY(yb.index()) = *yb;
          }
        for (auto&& xyb : indexed(hxy))
          {
            CountsXY(xyb.index(0) * nbins + xyb.index(1)) = *xyb;
          }

        double H_X = shan_entropy(CountsX);
        double H_Y = shan_entropy(CountsY);
        double H_XY = shan_entropy(CountsXY);
        //std::cout << "Hist: "<< H_X << " " << H_Y << " "<< H_XY << std::endl;
        //std::ofstream outfile("hist.out");
        //std::copy(xycount.begin(),xycount.end(),std::ostream_iterator<double>(outfile,"\n"));

        double I = 2.0 * H_XY - H_X - H_Y;
        double H = H_X + H_Y - H_XY;
        //return std::exp(-H);
        return std::sqrt(I);
#endif
      }

    void GaussHist(const jif3D::rvec &x, const jif3D::rvec &y, double xmin, double xmax,
        double ymin, double ymax, size_t nbins, jif3D::rvec &CountsX,
        jif3D::rvec &CountsY, jif3D::rvec &CountsXY)
      {

        boost::math::normal norm;
        const size_t nparm = x.size();

        double xw = (xmax - xmin) / nbins;
        double yw = (ymax - ymin) / nbins;
        const double binwidth = std::min(xw, yw) / 2.0;

        CountsX.resize(nbins);
        CountsY.resize(nbins);
        CountsXY.resize(nbins * nbins);
        std::fill(CountsX.begin(), CountsX.end(), 0.0);
        std::fill(CountsY.begin(), CountsY.end(), 0.0);
        std::fill(CountsXY.begin(), CountsXY.end(), 0.0);
        //std::ofstream outfile("data.out");

        for (size_t i = 0; i < nparm; ++i)
          {
            //outfile << x(i) << " " << y(i) << std::endl;
            //int xc = std::floor(std::min((x(i) - xmin) / xw, nbins - 1));
            //int yc = std::floor(std::min((y(i) - ymin) / yw, nbins - 1));
            for (size_t j = 0; j < nbins; ++j)
              {
                const double distx = x(i) - (xmin + (j + 0.5) * xw);
                for (size_t k = 0; k < nbins; ++k)
                  {
                    const double disty = y(i) - (ymin + (k + 0.5) * yw);
                    const double val = boost::math::pdf(norm,
                        std::sqrt(distx * distx + disty * disty) / binwidth);
                    CountsXY(j * nbins + k) += val;
                    CountsX(j) += val;
                    CountsY(k) += val;
                  }
              }
          }
      }

    jif3D::rvec diff_MIGauss(const jif3D::rvec &x, const jif3D::rvec &y, double xmin,
        double xmax, double ymin, double ymax, size_t nbins, const jif3D::rvec &CountsX,
        const jif3D::rvec &CountsY, const jif3D::rvec &CountsXY)
      {
        boost::math::normal norm;
        const size_t nparm = x.size();
        double xw = (xmax - xmin) / nbins;
        double yw = (ymax - ymin) / nbins;
        const double binwidth = std::min(xw, yw) / 2.0;

        jif3D::rvec Result(2 * nparm, 0.0);
        jif3D::rvec dsx = diff_shan_entropy(CountsX);
        jif3D::rvec dsy = diff_shan_entropy(CountsY);
        jif3D::rvec dsxy = diff_shan_entropy(CountsXY);
        for (size_t i = 0; i < nparm; ++i)
          {
            //int xc = std::floor(std::min((x(i) - xmin) / xw, nbins - 1));
            //int yc = std::floor(std::min((y(i) - ymin) / yw, nbins - 1));
            for (size_t j = 0; j < nbins; ++j)
              {
                const double distx = x(i) - (xmin + (j + 0.5) * xw);
                for (size_t k = 0; k < nbins; ++k)
                  {
                    const double disty = y(i) - (ymin + (k + 0.5) * yw);
                    const double dist = std::sqrt(distx * distx + disty * disty);
                    const double val = boost::math::pdf(norm, dist / binwidth)
                        / (binwidth * binwidth);
                    Result(i) -= 2.0 * distx * val * dsxy(j * nbins + k);
                    Result(i + nparm) -= 2.0 * disty * val * dsxy(j * nbins + k);
                    Result(i) += distx * val * dsx(j);
                    Result(i + nparm) += disty * val * dsx(j);
                    Result(i) += distx * val * dsy(k);
                    Result(i + nparm) += disty * val * dsy(k);
                  }
              }
          }

        return Result;
      }

    double MIGauss(const jif3D::rvec &x, const jif3D::rvec &y, double xmin, double xmax,
        double ymin, double ymax, size_t nbins, jif3D::rvec &CountsX,
        jif3D::rvec &CountsY, jif3D::rvec &CountsXY)
      {
        CountsXY.resize(nbins * nbins);
        CountsX.resize(nbins);
        CountsY.resize(nbins);
        GaussHist(x, y, xmin, xmax, ymin, ymax, nbins, CountsX, CountsY, CountsXY);
        double H_X = shan_entropy(CountsX);
        double H_Y = shan_entropy(CountsY);
        double H_XY = shan_entropy(CountsXY);
        //std::cout << "Gauss: " << H_X << " " << H_Y << " " << H_XY << std::endl;
        //std::ofstream outfile("gauss.out");
        //std::copy(CountsXY.begin(), CountsXY.end(),
        //    std::ostream_iterator<double>(outfile, "\n"));
        double I = 2.0 * H_XY - H_X - H_Y;
        double H = H_X + H_Y - H_XY;

        return std::sqrt(I);

      }

//! The implementation of the cross-gradient calculation
    void MutualInformationConstraint::ImplDataDifference(const jif3D::rvec &Model,
        jif3D::rvec &Diff)
      {

        Diff.resize(1);
        const size_t nparm = Model.size() / 2;
        jif3D::rvec x = ublas::subrange(Model, 0, nparm);
        jif3D::rvec y = ublas::subrange(Model, nparm, 2 * nparm);
        double H = 0;
        if (Calc == gauss)
          {
            H = MIGauss(x, y, xmin, xmax, ymin, ymax, nbins, CountsX, CountsY, CountsXY);
          }
        else
          {
            H = MIHist(x, y, xmin, xmax, ymin, ymax, nbins, CountsX, CountsY, CountsXY);
          }
        std::string name = "mi" + std::to_string(++eval) + ".nc";
        plothist(name, xmin, xmax, ymin, ymax, nbins, CountsXY);
        //double H2 = MIHist(x, y, xmin, xmax, ymin, ymax, nbins);
        //std::cout << "Gauss: " << H << " Hist: " << H2 << std::endl;
        Diff(0) = H;
      }

//! The gradient of the cross-gradient objective function with respect to the model parameters
    jif3D::rvec MutualInformationConstraint::ImplGradient(const jif3D::rvec &Model,
        const jif3D::rvec &Diff)
      {
        const size_t nparm = Model.size() / 2;
        jif3D::rvec x = ublas::subrange(Model, 0, nparm);
        jif3D::rvec y = ublas::subrange(Model, nparm, 2 * nparm);

        const size_t nmodel = Model.size();
        jif3D::rvec result = diff_MIGauss(x, y, xmin, xmax, ymin, ymax, nbins, CountsX,
            CountsY, CountsXY);
        return result;
      }

    MutualInformationConstraint::MutualInformationConstraint(double min1, double max1,
        double min2, double max2, size_t nb, MICalcType C) :
        CountsXY(), CountsX(), CountsY(), xmin(min1), xmax(max1), ymin(min2), ymax(max2), nbins(
            nb), eval(0), Calc(C)
      {
        // TODO Auto-generated constructor stub
        //mine = boost::make_shared<MINE>(0.6, 15, EST_MIC_APPROX);
      }

    MutualInformationConstraint::~MutualInformationConstraint()
      {
        // TODO Auto-generated destructor stub
      }

  } /* namespace jif3D */
