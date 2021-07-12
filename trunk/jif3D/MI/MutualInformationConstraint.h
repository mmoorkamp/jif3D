/*
 * MutualInformationConstraint.h
 *
 *  Created on: 2 Jul 2019
 *      Author: moorkamp
 */

#ifndef MI_MUTUALINFORMATIONCONSTRAINT_H_
#define MI_MUTUALINFORMATIONCONSTRAINT_H_

#include "../Global/Serialization.h"
#include "../Inversion/ObjectiveFunction.h"
#include <boost/shared_ptr.hpp>

namespace jif3D
  {

    double shan_entropy(const jif3D::rvec &x);
    jif3D::rvec diff_shan_entropy(const jif3D::rvec &x);

    class MutualInformationConstraint: public ObjectiveFunction
      {
    private:
      jif3D::rvec CountsXY;
      jif3D::rvec CountsX;
      jif3D::rvec CountsY;
      double xmin;
      double xmax;
      double ymin;
      double ymax;
      size_t nbins;
      double xw;
      double yw;
      double binwidth;
      double step;
      const size_t ngauss = 100000;
      jif3D::rvec GaussVals;
      inline double InterGauss(double val)
        {
          const int index = std::round(val / step);
          return index >= ngauss ? 0 : GaussVals(index);
        }
      double MIGauss(const jif3D::rvec &x, const jif3D::rvec &y, double xmin, double xmax,
          double ymin, double ymax, size_t nbins, jif3D::rvec &CountsX,
          jif3D::rvec &CountsY, jif3D::rvec &CountsXY);
      void GaussHist(const jif3D::rvec &x, const jif3D::rvec &y, double xmin, double xmax,
          double ymin, double ymax, size_t nbins, jif3D::rvec &CountsX,
          jif3D::rvec &CountsY, jif3D::rvec &CountsXY);
      jif3D::rvec diff_MIGauss(const jif3D::rvec &x, const jif3D::rvec &y, double xmin,
          double xmax, double ymin, double ymax, size_t nbins, const jif3D::rvec &CountsX,
          const jif3D::rvec &CountsY, const jif3D::rvec &CountsXY);
      //! The implementation of the MI calculation
      virtual void
      ImplDataDifference(const jif3D::rvec &Model, jif3D::rvec &Diff) override;
      //! The gradient of the cross-gradient objective function with respect to the model parameters
      virtual jif3D::rvec ImplGradient(const jif3D::rvec &Model, const jif3D::rvec &Diff)
          override;
    public:
      //! The clone function provides a virtual constructor
      virtual MutualInformationConstraint* clone() const override
        {
          return new MutualInformationConstraint(*this);
        }

      MutualInformationConstraint(double min1, double max1, double min2, double max2,
          size_t nb);
      virtual ~MutualInformationConstraint();
      };

  } /* namespace jif3D */

#endif /* MI_MUTUALINFORMATIONCONSTRAINT_H_ */
