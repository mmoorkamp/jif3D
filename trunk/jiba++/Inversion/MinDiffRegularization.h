//============================================================================
// Name        : MinDiffRegularization.h
// Author      : May 7, 2009
// Version     :
// Copyright   : 2009, mmoorkamp
//============================================================================


#ifndef MINDIFFREGULARIZATION_H_
#define MINDIFFREGULARIZATION_H_
#include "ObjectiveFunction.h"

namespace jiba
  {

    class MinDiffRegularization: public ObjectiveFunction
      {
    private:
      jiba::rvec Reference;
      virtual void ImplDataDifference(const jiba::rvec &Model, jiba::rvec &Diff)
        {
          Diff = Model - Reference;
        }
      //! The abstract interface for the gradient calculation
      virtual jiba::rvec ImplGradient(const jiba::rvec &Model,
          const jiba::rvec &Diff)
        {
          return 2.0*Diff;
        }
    public:
      void SetReferenceModel(const jiba::rvec &Model)
        {
          Reference = Model;
        }
      MinDiffRegularization();
      virtual ~MinDiffRegularization();
      };

  }

#endif /* MINDIFFREGULARIZATION_H_ */
