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
    public:
      //! The clone function provides a virtual constructor
      virtual MutualInformationConstraint *clone() const override
        {
          return new MutualInformationConstraint(*this);
        }
      //! The implementation of the cross-gradient calculation
      virtual void
      ImplDataDifference(const jif3D::rvec &Model, jif3D::rvec &Diff) override;
      //! The gradient of the cross-gradient objective function with respect to the model parameters
      virtual jif3D::rvec ImplGradient(const jif3D::rvec &Model, const jif3D::rvec &Diff)
          override;
      MutualInformationConstraint(double min1, double max1, double min2, double max2, size_t nb);
      virtual ~MutualInformationConstraint();
      };

  } /* namespace jif3D */

#endif /* MI_MUTUALINFORMATIONCONSTRAINT_H_ */
