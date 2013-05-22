//============================================================================
// Name        : LimitedMemoryQuasiNewton.h
// Author      : Apr 20, 2009
// Version     :
// Copyright   : 2009, mmoorkamp
//============================================================================


#ifndef LIMITEDMEMORYQUASINEWTON_H_
#define LIMITEDMEMORYQUASINEWTON_H_

#include "GradientBasedOptimization.h"
#include <boost/shared_ptr.hpp>
#include <vector>

namespace jif3D
  {
    /** \addtogroup inversion General routines for inversion */
    /* @{ */
    //! An implementation of the L-BFGS limited memory quasi-Newton optimization technique
    /*! This class implements the L-BFGS optimization method as described in Nocedal and Wright,
     * "Numerical Optimization", p. 225 with additions to consider model covariance as described
     * in Tarantolla for the quasi-Newton approach.
     *
     * We have to specify a number of correction pairs when constructing an object. For each correction pair
     * we need 2*M additional storage, where M is the number of model parameters.
     */
    class LimitedMemoryQuasiNewton: public jif3D::GradientBasedOptimization
      {
    private:
      //! the stepsize
      double mu;
      //! The maximum number of iterations in the line search
      size_t LineIter;
      //! The maximum number of correction pairs
      const size_t MaxPairs;
      //! The difference between the last MaxPairs updates, \f$s_k\f$ in eq. 9.4 of Nocedal and Wright
      std::vector<boost::shared_ptr<jif3D::rvec> > SHistory;
      //! The difference between the last MaxPairs gradients, \f$y_k\f$ in eq. 9.4 of Nocedal and Wright
      std::vector<boost::shared_ptr<jif3D::rvec> > YHistory;
      //! The implementation of a single step
      virtual void StepImplementation(jif3D::rvec &CurrentModel);
    public:
      //! Set the maximum number of forward evaluations the inversion spends in one line-search
      void SetLineSearchIterations(const size_t iter)
        {
          LineIter = iter;
        }
      //! The constructor needs an objective function object and the maximum number of correction pairs
      /*! When constructing the object we need two important pieces of information, the objective function
       * we want to minimize and the maximum number of correction pairs.
       * @param ObjFunction A shared pointer to an objective function object
       * @param n The maximum number of correction pairs
       */
      explicit LimitedMemoryQuasiNewton(boost::shared_ptr<
          jif3D::ObjectiveFunction> ObjFunction, const size_t n = 5);
      virtual ~LimitedMemoryQuasiNewton();
      };
  /* @} */
  }

#endif /* LIMITEDMEMORYQUASINEWTON_H_ */
