//============================================================================
// Name        : NonLinearConjugateGradient.h
// Author      : Apr 20, 2009
// Version     :
// Copyright   : 2009, mmoorkamp
//============================================================================


#ifndef NONLINEARCONJUGATEGRADIENT_H_
#define NONLINEARCONJUGATEGRADIENT_H_

#include "NonLinearOptimization.h"

namespace jiba
  {

    class NonLinearConjugateGradient: public jiba::NonLinearOptimization
      {
      private:
        jiba::rvec OldGradient;
        jiba::rvec OldDirection;
        double OldOmega;
        double mu;
        virtual void StepImplementation(jiba::rvec &CurrentModel);
    public:
      NonLinearConjugateGradient(boost::shared_ptr<jiba::ObjectiveFunction> ObjFunction);
      virtual ~NonLinearConjugateGradient();
      };

  }

#endif /* NONLINEARCONJUGATEGRADIENT_H_ */
