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
#include "../Global/Jif3DGlobal.h"

namespace jif3D
  {

    class J3DEXPORT MultiplicativeObjective : public JointObjective
      {
      private:
      //the implementation of the misfit calculation
      virtual void
      ImplDataDifference(const jif3D::rvec &Model, jif3D::rvec &Diff);
      //the implementation of the gradient calculation
      virtual jif3D::rvec ImplGradient(const jif3D::rvec &Model, const jif3D::rvec &Diff);
    public:
      MultiplicativeObjective(bool Verbose = false);
      virtual ~MultiplicativeObjective();
      };

  } /* namespace jif3D */
#endif /* MULTIPLICATIVEOBJECTIVE_H_ */
