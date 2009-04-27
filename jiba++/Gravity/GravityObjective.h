//============================================================================
// Name        : GravityObjective.h
// Author      : Apr 16, 2009
// Version     :
// Copyright   : 2009, mmoorkamp
//============================================================================


#ifndef GRAVITYOBJECTIVE_H_
#define GRAVITYOBJECTIVE_H_

#include "../Inversion/ObjectiveFunction.h"
#include "ThreeDGravityCalculator.h"
#include <boost/shared_ptr.hpp>

namespace jiba
  {
    /** \addtogroup gravity Gravity forward modeling, display and inversion */
    /* @{ */
    class GravityObjective: public ObjectiveFunction
      {
    private:
      boost::shared_ptr<jiba::ThreeDGravityCalculator> Calculator;
      jiba::rvec ObservedData;
      jiba::ThreeDGravityModel DensityModel;
      virtual void ImplDataDifference(const jiba::rvec &Model, jiba::rvec &Diff);
      virtual jiba::rvec ImplGradient(const jiba::rvec &Model, const  jiba::rvec &Diff);
      virtual void SetDataTransformAction()
        {
          Calculator->SetDataTransform(GetDataTransform());
        }
    public:
      void SetObservedData(const jiba::rvec &Data)
        {
          ObservedData = Data;
        }
      void SetModelGeometry(const jiba::ThreeDGravityModel &Model)
        {
          DensityModel = Model;
        }
      GravityObjective(bool ftg = false, bool cuda = false);
      virtual ~GravityObjective();
      };
  /* @} */
  }

#endif /* GRAVITYOBJECTIVE_H_ */
