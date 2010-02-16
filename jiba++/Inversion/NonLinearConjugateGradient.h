//============================================================================
// Name        : NonLinearConjugateGradient.h
// Author      : Apr 20, 2009
// Version     :
// Copyright   : 2009, mmoorkamp
//============================================================================


#ifndef NONLINEARCONJUGATEGRADIENT_H_
#define NONLINEARCONJUGATEGRADIENT_H_

#include "GradientBasedOptimization.h"

namespace jiba
  {
    /** \addtogroup inversion General routines for inversion */
    /* @{ */
    //! An implementation of the non-linear conjugate gradient (NLCG) algorithm
    /*! This class implements the NLCG algorithm as desribed in Tarantola.
     * It only requires the additional storage of the last model update
     * and the last gradient. Some of our tests indicate however that L-BFGS
     * performs better at a moderate cost of additional storage.
     */
    class NonLinearConjugateGradient: public jiba::GradientBasedOptimization
      {
    private:
      //! The gradient at the last iteration
      jiba::rvec OldGradient;
      //! The last search direction
      jiba::rvec OldDirection;
      //! See Tarantolla
      double OldOmega;
      //! The step-size
      double mu;
      //! Implements a single NLCG step including line search
      virtual void StepImplementation(jiba::rvec &CurrentModel);
    public:
      //! We only have to provide a shared pointer to an objective function object
      /*! The constructor needs to know what to minimize
       * @param ObjFunction A shared pointer to an objective function obejct
       */
      explicit NonLinearConjugateGradient(boost::shared_ptr<
          jiba::ObjectiveFunction> ObjFunction);
      virtual ~NonLinearConjugateGradient();
      };
  /* @} */
  }

#endif /* NONLINEARCONJUGATEGRADIENT_H_ */
