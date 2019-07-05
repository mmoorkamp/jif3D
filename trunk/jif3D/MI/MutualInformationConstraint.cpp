/*
 * MutualInformationConstraint.cpp
 *
 *  Created on: 2 Jul 2019
 *      Author: moorkamp
 */

#include "MutualInformationConstraint.h"
#include <omp.h>
#include <boost/histogram.hpp>
#include <boost/algorithm/minmax_element.hpp>
#include <algorithm>

namespace jif3D
  {

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
            if (x(i) > 0.0)
              {
                const double val = x(i) / sum;
                H -= val * std::log(val);
              }
          }
        return H;
      }

    double CalcMI(const jif3D::rvec &x, const jif3D::rvec &y, double xmin, double xmax,
        double ymin, double ymax, size_t nbins)
      {
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
        jif3D::rvec xcount(nbins), ycount(nbins), xycount(nbins * nbins);
        for (auto xb : indexed(hx))
          {
            xcount(xb.index()) = *xb;
          }
        for (auto yb : indexed(hy))
          {
            ycount(yb.index()) = *yb;
          }
        for (auto xyb : indexed(hxy))
          {
            xycount(xyb.index(0) * nbins + xyb.index(1)) = *xyb;
          }

        double H_X = shan_entropy(xcount);
        double H_Y = shan_entropy(ycount);
        double H_XY = shan_entropy(xycount);
        //std::cout << H_X << " " << H_Y << " "<< H_XY << std::endl;
        //double H = 2.0 * H_XY - H_X - H_Y;
        double H = H_X + H_Y - H_XY;
        return std::exp(-H);
      }

    //! The implementation of the cross-gradient calculation
    void MutualInformationConstraint::ImplDataDifference(const jif3D::rvec &Model,
        jif3D::rvec &Diff)
      {
        // enables _c suffix
        Diff.resize(1);
        const size_t nparm = Model.size() / 2;
        jif3D::rvec x = ublas::subrange(Model, 0, nparm);
        jif3D::rvec y = ublas::subrange(Model, nparm, 2 * nparm);
//        auto xmm = boost::minmax_element(x.begin(), x.end());
//        auto ymm = boost::minmax_element(y.begin(), y.end());
//        double eps = std::max(*(xmm.second), *(ymm.second)) * 0.002;
//        double xmin = *(xmm.first) + eps;
//        double xmax = *(xmm.second) - eps;
//        double ymin = *(ymm.first) + eps;
//        double ymax = *(ymm.second) - eps;
        const size_t nbins = 100;
        double H = CalcMI(x, y, xmin, xmax, ymin, ymax, nbins);
        Diff(0) = H;
      }

    //! The gradient of the cross-gradient objective function with respect to the model parameters
    jif3D::rvec MutualInformationConstraint::ImplGradient(const jif3D::rvec &Model,
        const jif3D::rvec &Diff)
      {
        const size_t nparm = Model.size() / 2;
        jif3D::rvec x = ublas::subrange(Model, 0, nparm);
        jif3D::rvec y = ublas::subrange(Model, nparm, 2 * nparm);
//        auto xmm = boost::minmax_element(x.begin(), x.end());
//        auto ymm = boost::minmax_element(y.begin(), y.end());
//        double eps = std::max(*(xmm.second), *(ymm.second)) * 0.002;
//        double xmin = *(xmm.first) + eps;
//        double xmax = *(xmm.second) - eps;
//        double ymin = *(ymm.first) + eps;
//        double ymax = *(ymm.second) - eps;
        const size_t nbins = 100;

        const size_t nmodel = Model.size();
        jif3D::rvec result(nmodel, 0.0);
        omp_lock_t lck;
        omp_init_lock(&lck);
//#pragma omp parallel for
        for (size_t i = 0; i < nmodel; ++i)
          {
            double delta = Model(i) * 0.001;

            jif3D::rvec TmpModel(Model);
            TmpModel(i) += delta;
            double H = CalcMI(ublas::subrange(TmpModel, 0, nparm),
                ublas::subrange(TmpModel, nparm, 2 * nparm), xmin, xmax, ymin, ymax,
                nbins);
            double Forward = jif3D::pow2(H);
            TmpModel(i) -= 2 * delta;
            H = CalcMI(ublas::subrange(TmpModel, 0, nparm),
                ublas::subrange(TmpModel, nparm, 2 * nparm), xmin, xmax, ymin, ymax,
                nbins);
            double Backward = jif3D::pow2(H);
            //omp_set_lock(&lck);
            result(i) = (Forward - Backward) / (2.0 * delta);
            //omp_unset_lock(&lck);
          }

        omp_destroy_lock(&lck);
        return result;
      }

    MutualInformationConstraint::MutualInformationConstraint(double min1, double max1,
        double min2, double max2) :
        xmin(min1), xmax(max1), ymin(min2), ymax(max2)
      {
        // TODO Auto-generated constructor stub
        //mine = boost::make_shared<MINE>(0.6, 15, EST_MIC_APPROX);
      }

    MutualInformationConstraint::~MutualInformationConstraint()
      {
        // TODO Auto-generated destructor stub
      }

  } /* namespace jif3D */
