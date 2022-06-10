/*
 * MutualInformationConstraint.cpp
 *
 *  Created on: 2 Jul 2019
 *      Author: moorkamp
 */

#include "MutualInformationConstraint.h"
#include "Entropy.h"
#include "../Global/NetCDFTools.h"
#include "../Global/NumUtil.h"
#include "../Global/FatalException.h"
#include <omp.h>

#include <boost/algorithm/minmax_element.hpp>
#include <boost/math/distributions/normal.hpp>
#include <boost/math/special_functions/fpclassify.hpp>
#include <algorithm>
#include <cmath>
#include <fstream>

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

    void MutualInformationConstraint::GaussHist(const jif3D::rvec &x,
        const jif3D::rvec &y, double xmin, double xmax, double ymin, double ymax,
        size_t nbins, jif3D::rvec &CountsX, jif3D::rvec &CountsY, jif3D::rvec &CountsXY)
      {

        boost::math::normal norm;
        const size_t nparm = x.size();

        CountsX.resize(nbins);
        CountsY.resize(nbins);
        CountsXY.resize(nbins * nbins);
        std::fill(CountsX.begin(), CountsX.end(), 0.0);
        std::fill(CountsY.begin(), CountsY.end(), 0.0);
        std::fill(CountsXY.begin(), CountsXY.end(), 0.0);

        double factor = 1.0 / nparm;
        //int xc = std::floor(std::min((x(i) - xmin) / xw, nbins - 1));
        //int yc = std::floor(std::min((y(i) - ymin) / yw, nbins - 1));
#pragma omp parallel for
        for (size_t j = 0; j < nbins; ++j)
          {
            for (size_t i = 0; i < nparm; ++i)
              {
                if (ValidityVector(i) >= 1.0)
                  {
                    const double distx = x(i) - (xmin + (j + 0.5) * xw);
                    for (size_t k = 0; k < nbins; ++k)
                      {
                        const double disty = y(i) - (ymin + (k + 0.5) * yw);
                        const double dist = std::sqrt(distx * distx + disty * disty)
                            / stddev;

                        const double val = boost::math::pdf(norm, dist) * factor;
                        //const double val = InterGauss(dist);

                        CountsXY(j * nbins + k) += val;
                        CountsX(j) += val;
#pragma omp atomic
                        CountsY(k) += val;
                      }
                  }
              }
          }
      }

    jif3D::rvec MutualInformationConstraint::diff_MIGauss(const jif3D::rvec &x,
        const jif3D::rvec &y, double xmin, double xmax, double ymin, double ymax,
        size_t nbins, const jif3D::rvec &CountsX, const jif3D::rvec &CountsY,
        const jif3D::rvec &CountsXY)
      {
        boost::math::normal norm;
        const size_t nparm = x.size();

        jif3D::rvec Result(2 * nparm, 0.0);
        jif3D::rvec dsx = diff_shan_entropy(CountsX);
        jif3D::rvec dsy = diff_shan_entropy(CountsY);
        jif3D::rvec dsxy = diff_shan_entropy(CountsXY);
        double factor = 1.0 / nparm;

#pragma omp parallel for
        for (size_t i = 0; i < nparm; ++i)
          {
            if (ValidityVector(i) >= 1.0)
              {
                //int xc = std::floor(std::min((x(i) - xmin) / xw, nbins - 1));
                //int yc = std::floor(std::min((y(i) - ymin) / yw, nbins - 1));
                for (size_t j = 0; j < nbins; ++j)
                  {
                    const double distx = x(i) - (xmin + (j + 0.5) * xw);
                    for (size_t k = 0; k < nbins; ++k)
                      {
                        const double disty = y(i) - (ymin + (k + 0.5) * yw);
                        const double dist = std::sqrt(distx * distx + disty * disty)
                            / stddev;
                        const double val = boost::math::pdf(norm, dist)
                            / (stddev * stddev) * factor;
                        //const double val = InterGauss(dist)/ (binwidth * binwidth);
                        Result(i) -= 2.0 * distx * val * dsxy(j * nbins + k);
                        Result(i + nparm) -= 2.0 * disty * val * dsxy(j * nbins + k);
                        Result(i) += distx * val * dsx(j);
                        Result(i + nparm) += disty * val * dsx(j);
                        Result(i) += distx * val * dsy(k);
                        Result(i + nparm) += disty * val * dsy(k);
                      }
                  }
              }
          }

        return Result;
      }

    double MutualInformationConstraint::MIGauss(const jif3D::rvec &x,
        const jif3D::rvec &y, double xmin, double xmax, double ymin, double ymax,
        size_t nbins, jif3D::rvec &CountsX, jif3D::rvec &CountsY, jif3D::rvec &CountsXY)
      {
        CountsXY.resize(nbins * nbins);
        CountsX.resize(nbins);
        CountsY.resize(nbins);
        GaussHist(x, y, xmin, xmax, ymin, ymax, nbins, CountsX, CountsY, CountsXY);
        double H_X = shan_entropy(CountsX);
        double H_Y = shan_entropy(CountsY);
        double H_XY = shan_entropy(CountsXY);

        double I = 2.0 * H_XY - H_X - H_Y;
        //double H = H_X + H_Y - H_XY;

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
        if (ValidityVector.empty())
          {
            ValidityVector.resize(nparm);
            std::fill(ValidityVector.begin(), ValidityVector.end(), 1.0);
          }
        else
          {
            if (ValidityVector.size() != nparm)
              {
                throw jif3D::FatalException(
                    "Validity vector size: " + std::to_string(ValidityVector.size())
                        + " does not match grid: " + std::to_string(nparm));
              }
          }
        double H = 0;

        H = MIGauss(x, y, xmin, xmax, ymin, ymax, nbins, CountsX, CountsY, CountsXY);
        //std::string name = "mi.nc";
        //plothist(name, xmin, xmax, ymin, ymax, nbins, CountsXY);
        //std::cout << "Min x: " << *std::min_element(x.begin(),x.end()) << std::endl;
        //std::cout << "Max x: " << *std::max_element(x.begin(),x.end()) << std::endl;
        //std::cout << "Min y: " << *std::min_element(y.begin(),y.end()) << std::endl;
        //std::cout << "Max y: " << *std::max_element(y.begin(),y.end()) << std::endl;
        Diff(0) = H;
      }

//! The gradient of the cross-gradient objective function with respect to the model parameters
    jif3D::rvec MutualInformationConstraint::ImplGradient(const jif3D::rvec &Model,
        const jif3D::rvec &Diff)
      {
        const size_t nparm = Model.size() / 2;
        jif3D::rvec x = ublas::subrange(Model, 0, nparm);
        jif3D::rvec y = ublas::subrange(Model, nparm, 2 * nparm);
        jif3D::rvec result = diff_MIGauss(x, y, xmin, xmax, ymin, ymax, nbins, CountsX,
            CountsY, CountsXY);
        return result;
      }

    MutualInformationConstraint::MutualInformationConstraint(double min1, double max1,
        double min2, double max2, size_t nb, const jif3D::rvec &ValVec) :
        CountsXY(), CountsX(), CountsY(), ValidityVector(ValVec), xmin(min1), xmax(max1), ymin(min2), ymax(max2), nbins(
            nb)
      {
        boost::math::normal norm;
        const double maxdist = std::sqrt(
            jif3D::pow2(xmax - xmin) + jif3D::pow2(ymax - ymin));
        xw = (xmax - xmin) / nbins;
        yw = (ymax - ymin) / nbins;
        stddev = std::min(xw, yw) / 2.0;
        // TODO Auto-generated constructor stub
        //mine = boost::make_shared<MINE>(0.6, 15, EST_MIC_APPROX);
      }

    MutualInformationConstraint::~MutualInformationConstraint()
      {
        // TODO Auto-generated destructor stub
      }

  } /* namespace jif3D */
