//============================================================================
// Name        : NonLinearOptimization.h
// Author      : Apr 16, 2009
// Version     :
// Copyright   : 2009, mmoorkamp
//============================================================================

#ifndef NONLINEAROPTIMIZATION_H_
#define NONLINEAROPTIMIZATION_H_
#ifdef HAVEHPX
#include <hpx/config.hpp>
#endif
#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>
#include "ObjectiveFunction.h"
#include "GeneralCovariance.h"
#include "DiagonalCovariance.h"
#include "../Global/VecMat.h"
#include "../Global/Jif3DGlobal.h"

namespace jif3D
  {
    /** \addtogroup inversion General routines for inversion */
    /* @{ */
    //! This is the base class for all non-linear optimization methods.
    /*! This class provides a least common denominator interface for
     * all non-linear optimization methods. It also stores some basic
     * information about covariances and misfit that can be used in all
     * types of optimization methods. Derived classes habe to implement
     * two virtual functions, EvaluateModel and StepImplementation. See
     * the respective description for requirements and responsibilities.
     *
     * The main method MakeStep performs a single optimization step.
     *
     */
    class J3DEXPORT NonLinearOptimization
      {
    private:
      //! Evaluate all necessary information for the current model
      /*! The implementation of this method has to provide all required
       * information about the current model for the optimization algorithm.
       * For example the misift and the gradient o0f the current model for gradient based methods
       * \see GradientBasedOptimization.
       * @param CurrentModel The vector of current model parameters
       */
      virtual void EvaluateModel(const jif3D::rvec &CurrentModel) = 0;
      //! Perform a single step of the optimization
      /*! This is the main implementation of the optimization algorithm. It needs to update
       * the current model vector by performing a single step in the optimization. The implementation
       * can rely on the fact that EvaluateModel has been called before this routine and the misfit
       * information for the current model and everything else that is updated there is correct.
       * @param CurrentModel The current model, contains the updated model on exit
       */
      virtual void StepImplementation(jif3D::rvec &CurrentModel) = 0;
      //! The object applying the model covariance
      boost::shared_ptr<jif3D::GeneralCovariance> Covar;
      //! The objective function object
      boost::shared_ptr<jif3D::ObjectiveFunction> Objective;
    protected:
      //! The current misfit of the objective function
      double Misfit;
      //! Direct access to the objective function object for derived classes
      jif3D::ObjectiveFunction &GetObjective()
        {
          return *Objective.get();
        }
    public:
      const boost::shared_ptr<jif3D::GeneralCovariance> GetCovObj() const
        {
          return Covar;
        }
      void SetCovObj(boost::shared_ptr<jif3D::GeneralCovariance> Cv)
        {
          Covar = Cv;
        }
      //! Get the misfit after calling MakeStep
      double GetMisfit()
        {
          return Misfit;
        }
      //! Update the current model by making a single step with the optimization method
      /*! This function provides the outline for making a non-linear optimization step, evaluating the
       * model and calling the implementation of the step.
       * @param CurrentModel The current model vector, contains the updated model on exit.
       */
      void MakeStep(jif3D::rvec &CurrentModel);
      //! The constructor needs the objective function object, without it optimization does not make much sense
      explicit NonLinearOptimization(
          boost::shared_ptr<jif3D::ObjectiveFunction> ObjFunction,
          boost::shared_ptr<jif3D::GeneralCovariance> Cv = boost::make_shared<
              jif3D::DiagonalCovariance>());
      virtual ~NonLinearOptimization();
      };
  /* @} */
  }

#endif /* NONLINEAROPTIMIZATION_H_ */
