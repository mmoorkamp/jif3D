//============================================================================
// Name        : MultiplicativeObjective.h
// Author      : 21 Mar 2012
// Version     : 
// Copyright   : 2012, mm489
//============================================================================

#ifndef MULTIPLICATIVEOBJECTIVE_H_
#define MULTIPLICATIVEOBJECTIVE_H_

#include "JointObjective.h"
#include "ModelDistributor.h"
#include "../Global/FatalException.h"
#include <boost/shared_ptr.hpp>

namespace jiba
  {

    class MultiplicativeObjective : public JointObjective
      {
      private:
      //the implementation of the misfit calculation
      virtual void
      ImplDataDifference(const jiba::rvec &Model, jiba::rvec &Diff);
      //the implementation of the gradient calculation
      virtual jiba::rvec ImplGradient(const jiba::rvec &Model, const jiba::rvec &Diff);
    public:
      MultiplicativeObjective(bool Verbose = false);
      virtual ~MultiplicativeObjective();
      };

  } /* namespace jiba */
#endif /* MULTIPLICATIVEOBJECTIVE_H_ */
