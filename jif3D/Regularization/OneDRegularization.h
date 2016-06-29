//============================================================================
// Name        : OneDRegularization.h
// Author      : 14 Jun 2012
// Version     :
// Copyright   : 2012, mm489
//============================================================================

#ifndef ONEDREGULARIZATION_H_
#define ONEDREGULARIZATION_H_


#include "../Global/Serialization.h"
#include "../Global/VecMat.h"
#include "../Inversion/ObjectiveFunction.h"
#include "../Global/Jif3DGlobal.h"

namespace jif3D
  {
    /** \addtogroup Regularization classes to regularize the inversion */
    /* @{ */
    //! A first derivative based regularization method for 1D inversions
    /*! For our 1D experiments we need a special regularization class
     *
     */
    class J3DEXPORT OneDRegularization: public jif3D::ObjectiveFunction
      {
    private:
      //! The matrix we apply to the model vector to calculate the regularization
      comp_mat OperatorMatrix;
      //! The reference model substracted before calculating the regularization
      jif3D::rvec RefMod;
      friend class access;
      //! Provide serialization to be able to store objects and, more importantly for simpler MPI parallelization
      template<class Archive>
      void serialize(Archive & ar, const unsigned int version)
        {
          ar & base_object<ObjectiveFunction>(*this);
          ar & OperatorMatrix;
          ar & RefMod;
        }
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
      void SetReferenceModel(const jif3D::rvec &Model)
        {
          RefMod = Model;
        }
      virtual void ImplDataDifference(const jif3D::rvec &Model, jif3D::rvec &Diff);
      virtual jif3D::rvec ImplGradient(const jif3D::rvec &Model, const jif3D::rvec &Diff);
      OneDRegularization(const size_t nlayers);
      virtual ~OneDRegularization();
      };

  } /* namespace jif3D */
#endif /* ONEDREGULARIZATION_H_ */
