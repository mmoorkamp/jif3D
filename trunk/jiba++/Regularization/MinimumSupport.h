//============================================================================
// Name        : MinimumSupport.h
// Author      : Jul 9, 2012
// Version     :
// Copyright   : 2012, mmoorkamp
//============================================================================

#ifndef MINIMUMSUPPORT_H_
#define MINIMUMSUPPORT_H_

#include <boost/shared_ptr.hpp>
#include "MatOpRegularization.h"

namespace jiba
  {

    class MinimumSupport: public jiba::ObjectiveFunction
      {
    private:
      jiba::rvec RegDiff;
      double beta;
      boost::shared_ptr<jiba::MatOpRegularization> RegFunc;
    public:
      virtual void ImplDataDifference(const jiba::rvec &Model, jiba::rvec &Diff);
      virtual jiba::rvec ImplGradient(const jiba::rvec &Model, const jiba::rvec &Diff);
      MinimumSupport(boost::shared_ptr<jiba::MatOpRegularization> RF, double b = 1.0);
      virtual ~MinimumSupport();
      };

  } /* namespace jiba */
#endif /* MINIMUMSUPPORT_H_ */
