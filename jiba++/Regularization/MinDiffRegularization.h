//============================================================================
// Name        : MinDiffRegularization.h
// Author      : May 7, 2009
// Version     :
// Copyright   : 2009, mmoorkamp
//============================================================================


#ifndef MINDIFFREGULARIZATION_H_
#define MINDIFFREGULARIZATION_H_
#include "../Inversion/ObjectiveFunction.h"

namespace jiba
  {
    /** \addtogroup inversion General routines for inversion */
    /* @{ */
    //! Regularize by calculating the difference to a reference model
    /*! This class can be used to regularize the inversion by
     * looking for a model that is as close as possible to a reference model.
     */
    class MinDiffRegularization: public ObjectiveFunction
      {
    private:
      //! The reference model that we want to regularize to
      jiba::rvec Reference;
      //! The misfit is simply the difference between the vectors
      virtual void ImplDataDifference(const jiba::rvec &Model, jiba::rvec &Diff)
        {
          assert(Model.size() == Reference.size());
          Diff = Model - Reference;
        }
      //! The gradient is particularly simple
      virtual jiba::rvec ImplGradient(const jiba::rvec &Model,
          const jiba::rvec &Diff)
        {
          return 2.0 * Diff;
        }
    public:
      //! When constructing the object we have to specify a reference
      /*! The constructor always needs a reference model vector. This
       * has to have the same size as the model vector in the inversion.
       * @param Model The reference model
       */
      explicit MinDiffRegularization(const jiba::rvec &Model) :
        Reference(Model)
        {

        }
      virtual ~MinDiffRegularization()
        {

        }
      };
  /* @} */
  }

#endif /* MINDIFFREGULARIZATION_H_ */
