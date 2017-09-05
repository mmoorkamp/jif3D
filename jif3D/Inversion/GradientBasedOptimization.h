//============================================================================
// Name        : GradientBasedOptimization.h
// Author      : Jun 29, 2009
// Version     :
// Copyright   : 2009, mmoorkamp
//============================================================================


#ifndef GRADIENTBASEDOPTIMIZATION_H_
#define GRADIENTBASEDOPTIMIZATION_H_

#include "../Global/Jif3DGlobal.h"
#include "NonLinearOptimization.h"

namespace jif3D
  {
    /** \addtogroup inversion General routines for inversion */
    /* @{ */
    //! This class implements some functionality that is common to all gradient based optimization techniques
    /*! This class implements the virtual function EvaluateModel and provides storage for some quantities
     * that we always need in gradient based optimization.
     *
     * Important note: Typical line search algorithms have to evaluate the misfit and gradient at a number of trial
     * points including the model for the next iteration. In order to avoid duplicate misift and
     * gradient calculations we keep track of whether we have already calculated this information.
     * As we usually do not change the current model outside the inversion process we do this by
     * a simple boolean flag. This means that if you change any value in the current model
     * between two steps of the inversion you have to call InvalidateCache in order to ensure
     * proper calculation of misfit and gradient.
     */
    class J3DEXPORT GradientBasedOptimization : public jif3D::NonLinearOptimization
      {
    private:
      //! Have we evaluated the misfit and gradient for the current model?
      bool HaveEvaluated_;
      //! Calculate the misfit and the gradient for the current model
      virtual void EvaluateModel(const jif3D::rvec &CurrentModel) override;
    protected:
      //! Derived classes can call this function to indicate that they have updated the misfit and raw gradient
      void HaveEvaluated()
        {
          HaveEvaluated_ = true;
        }
      //! the gradient gamma  for the current model
      jif3D::rvec RawGrad;
      //! the direction of steepest ascent C_M gamma for the current model
      jif3D::rvec CovGrad;
      //! the search direction for the current opimization step
      jif3D::rvec SearchDir;
    public:
      //! Signal that the model has been changed and new misfit and gradient have to be calculated
      void InvalidateCache()
        {
          HaveEvaluated_ = false;
        }
      //! Get the norm of the current search direction
      double GetGradNorm()
        {
          return ublas::norm_2(SearchDir);
        }
      //! The constructor needs the objective function object, without it optimization does not make much sense
      explicit GradientBasedOptimization(boost::shared_ptr<
          jif3D::ObjectiveFunction> ObjFunction);
      virtual ~GradientBasedOptimization();
      };
  /* @} */
  }

#endif /* GRADIENTBASEDOPTIMIZATION_H_ */
