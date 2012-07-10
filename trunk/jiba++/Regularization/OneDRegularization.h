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
      comp_mat OperatorMatrix;
      jiba::rvec RefMod;
    public:
      //! Set the reference model for the roughness calculation, this is optional
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
