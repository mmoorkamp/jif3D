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
     * looking for a model that is as close as possible to a reference model. If no reference
     * model is specified it is assumed to be zero and therefore the smallest model
     * is sought.
     */
    class MinDiffRegularization: public ObjectiveFunction
      {
    private:
      //! The reference model that we want to regularize to
      jiba::rvec Reference;
      //! The misfit is simply the difference between the vectors
      virtual void ImplDataDifference(const jiba::rvec &Model, jiba::rvec &Diff)
        {
          if (Reference.empty())
            {
              Reference.resize(Model.size());
              std::fill(Reference.begin(), Reference.end(), 0.0);
            }
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
      //! Set the reference model for the roughness calculation, this is optional
      void SetReferenceModel(const jiba::rvec &Model)
        {
          Reference = Model;
        }
      MinDiffRegularization()
        {
        }
      virtual ~MinDiffRegularization()
        {
        }
      };
  /* @} */
  }

#endif /* MINDIFFREGULARIZATION_H_ */
