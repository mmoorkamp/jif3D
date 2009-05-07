//============================================================================
// Name        : TomographyObjective.h
// Author      : May 5, 2009
// Version     :
// Copyright   : 2009, mmoorkamp
//============================================================================


#ifndef TOMOGRAPHYOBJECTIVE_H_
#define TOMOGRAPHYOBJECTIVE_H_

#include "../Inversion/ObjectiveFunction.h"
#include "TomographyCalculator.h"
#include "ThreeDSeismicModel.h"

namespace jiba
  {

    class TomographyObjective: public ObjectiveFunction
      {
    private:
      jiba::ThreeDSeismicModel SlownessModel;
      jiba::rvec ObservedData;
      TomographyCalculator Calculator;
      virtual void
          ImplDataDifference(const jiba::rvec &Model, jiba::rvec &Diff);
      virtual jiba::rvec ImplGradient(const jiba::rvec &Model,
          const jiba::rvec &Diff);
      virtual void SetDataTransformAction()
        {
          //TODO Implement transformation if necessary
        }
    public:
      void SetObservedData(const jiba::rvec &Data)
        {
          ObservedData = Data;
        }
      void SetModelGeometry(const jiba::ThreeDSeismicModel &Model)
        {
          SlownessModel = Model;
        }
      TomographyObjective();
      virtual ~TomographyObjective();
      };

  }

#endif /* TOMOGRAPHYOBJECTIVE_H_ */
