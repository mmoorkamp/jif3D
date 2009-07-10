//============================================================================
// Name        : X3DObjective.h
// Author      : Jul 10, 2009
// Version     : 
// Copyright   : 2009, mmoorkamp
//============================================================================


#ifndef X3DOBJECTIVE_H_
#define X3DOBJECTIVE_H_

#include "../Inversion/ObjectiveFunction.h"
#include "X3DModel.h"

namespace jiba
  {

    class X3DObjective: public ObjectiveFunction
      {
    private:
      X3DModel ConductivityModel;
      jiba::rvec ObservedData;
      virtual void
      ImplDataDifference(const jiba::rvec &Model, jiba::rvec &Diff);
      //The implementation of the gradient calculation
      virtual jiba::rvec ImplGradient(const jiba::rvec &Model,
          const jiba::rvec &Diff);
      //! So far transformations have no effect for MT data !
      virtual void SetDataTransformAction()
        {
          //TODO Implement transformation if necessary
        }
    public:
      void SetObservedData(const jiba::rvec &Data)
        {
          ObservedData = Data;
        }
      void SetModelGeometry(const jiba::X3DModel &Model)
        {
          ConductivityModel = Model;
        }
      X3DObjective();
      virtual ~X3DObjective();
      };

  }

#endif /* X3DOBJECTIVE_H_ */
