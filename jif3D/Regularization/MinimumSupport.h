//============================================================================
// Name        : MinimumSupport.h
// Author      : Jul 9, 2012
// Version     :
// Copyright   : 2012, mmoorkamp
//============================================================================

#ifndef MINIMUMSUPPORT_H_
#define MINIMUMSUPPORT_H_

#include <boost/shared_ptr.hpp>
#include "RegularizationFunction.h"
#include "MatOpRegularization.h"

namespace jif3D
  {

    class MinimumSupport: public jif3D::RegularizationFunction
      {
    private:
      size_t ModelNumber;
      jif3D::rvec RegDiff;
      double beta;

      boost::shared_ptr<jif3D::MatOpRegularization> RegFunc;
      jif3D::ThreeDModelBase Geometry;
    public:
      //! The clone function provides a virtual constructor
      virtual MinimumSupport *clone() const
        {
          return new MinimumSupport(*this);
        }
      virtual void ImplDataDifference(const jif3D::rvec &Model, jif3D::rvec &Diff);
      virtual jif3D::rvec ImplGradient(const jif3D::rvec &Model, const jif3D::rvec &Diff);
      MinimumSupport(boost::shared_ptr<jif3D::MatOpRegularization> RF, double b = 1.0);
      virtual ~MinimumSupport();
      };

  } /* namespace jif3D */
#endif /* MINIMUMSUPPORT_H_ */
