//============================================================================
// Name        : OneDRegularization.h
// Author      : 14 Jun 2012
// Version     : 
// Copyright   : 2012, mm489
//============================================================================

#ifndef ONEDREGULARIZATION_H_
#define ONEDREGULARIZATION_H_

#include "../Global/VecMat.h"
#include "ObjectiveFunction.h"

namespace jiba
  {
    /** \addtogroup Regularization classes to regularize the inversion */
    /* @{ */
    //! A first derivative based regularization method for 1D inversions
    /*! For our 1D experiments we need a special regularization class
     *
     */
    class OneDRegularization: public jiba::ObjectiveFunction
      {
    private:
      //! The matrix we apply to the model vector to calculate the regularization
      comp_mat OperatorMatrix;
      //! The reference model substracted before calculating the regularization
      jiba::rvec RefMod;
    public:
      //! The clone function provides a virtual constructor
      virtual OneDRegularization *clone() const
        {
          return new OneDRegularization(*this);
        }
      //! Set the reference model for the roughness calculation, this is optional
      /*! We can optionally set a reference model. This will be substracted from
       * the current model before calculating the regularization value. If
       * the model is not set within this function it will be set to zero
       * automatically.
       * @param Model The model vector to be substracted from the current model.
       */
      void SetReferenceModel(const jiba::rvec &Model)
        {
          RefMod = Model;
        }
      virtual void ImplDataDifference(const jiba::rvec &Model, jiba::rvec &Diff);
      virtual jiba::rvec ImplGradient(const jiba::rvec &Model, const jiba::rvec &Diff);
      OneDRegularization(const size_t nlayers);
      virtual ~OneDRegularization();
      };

  } /* namespace jiba */
#endif /* ONEDREGULARIZATION_H_ */
