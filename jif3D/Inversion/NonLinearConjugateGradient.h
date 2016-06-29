//============================================================================
// Name        : NonLinearConjugateGradient.h
// Author      : Apr 20, 2009
// Version     :
// Copyright   : 2009, mmoorkamp
//============================================================================


#ifndef NONLINEARCONJUGATEGRADIENT_H_
#define NONLINEARCONJUGATEGRADIENT_H_

#include "GradientBasedOptimization.h"
#include "../Global/Jif3DGlobal.h"

namespace jif3D
  {
    /** \addtogroup inversion General routines for inversion */
    /* @{ */
    //! An implementation of the non-linear conjugate gradient (NLCG) algorithm
    /*! This class implements the NLCG algorithm as desribed in Tarantola.
     * It only requires the additional storage of the last model update
     * and the last gradient. Some of our tests indicate however that L-BFGS
     * performs better at a moderate cost of additional storage.
     */
    class J3DEXPORT NonLinearConjugateGradient: public jif3D::GradientBasedOptimization
      {
    private:
      //! The gradient at the last iteration
      jif3D::rvec OldGradient;
      //! The last search direction
      jif3D::rvec OldDirection;
      //! A scaling factor similar to L-BFGS to make help that after the first iteration we have the right step size
      double gamma;
      //! See Tarantolla
      double OldOmega;
      //! The step-size
      double mu;
      //! Implements a single NLCG step including line search
      virtual void StepImplementation(jif3D::rvec &CurrentModel);
    public:
      //! We only have to provide a shared pointer to an objective function object
      /*! The constructor needs to know what to minimize
       * @param ObjFunction A shared pointer to an objective function obejct
       */
      explicit NonLinearConjugateGradient(boost::shared_ptr<
          jif3D::ObjectiveFunction> ObjFunction);
      virtual ~NonLinearConjugateGradient();
      };
  /* @} */
  }

#endif /* NONLINEARCONJUGATEGRADIENT_H_ */
