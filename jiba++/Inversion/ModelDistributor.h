//============================================================================
// Name        : ModelDistributor.h
// Author      : Apr 15, 2009
// Version     :
// Copyright   : 2009, mmoorkamp
//============================================================================


#ifndef MODELDISTRIBUTOR_H_
#define MODELDISTRIBUTOR_H_

#include <cassert>
#include <vector>
#include <boost/function.hpp>
#include <boost/shared_ptr.hpp>
#include "../Global/VecMat.h"

namespace jiba
  {
    /** \addtogroup inversion General routines for inversion */
    /* @{ */

    //! This is the base class for all function objects that transform generalized model parameters to concrete model parameters
    /*! As described for the class ModelDistributor we need a flexible mechanism to decouple the inversion parameters from the
     * concrete physical parameters that each method works with. This base class provides the interface for all function objects
     * that implement such a transformation.
     */
    class GeneralModelTransform
      {
    public:
      GeneralModelTransform()
        {
        }
      virtual ~GeneralModelTransform()
        {
        }
      //! The declaration for the transformation
      /*! The transformation takes a generalized model vector and computes the physical model vector from it. This means
       * that the complete knowledge of the meaning of the generalized model parameters is encapsulated in this function.
       *
       * For example if the input parameter contains 3M elements that represent conductivity, density and velocity, the implementation
       * for the MT transformation only copies the conductivity parameters and returns a vector of length M, while the seismic transformation
       * takes the velocity parameters and transforms them into slowness.
       * @param FullModel The generalized model vector including all parameters
       * @return The physical quantities for one method
       */
      virtual jiba::rvec GeneralizedToPhysical(const jiba::rvec &FullModel) const = 0;
      //! Transform the physical model vector to a generalized model vector
      /*! This is the inverse transform to GeneralizedToPhysical. It is used
       * to transform a starting model that is given in terms of physical parameters
       * to generalized model parameters.
       * @param FullModel The physical model vector
       */
      virtual jiba::rvec PhysicalToGeneralized(const jiba::rvec &FullModel) const = 0;
      //! For the inversion we also need the derivative of the transformation
      /*! To calculate the derivative of the generalized parameters we need to apply the chain rule, as each method only gives
       * us the derivative with respect to the physical parameters. This function returns a vector with all the derivative values
       * for the current model.
       * @param FullModel The generalized model vector including all parameters
       * @param Derivative The gradient with respect to the untransformed parameters
       * @return The derivatives of the generalized parameters with respect to the physical parameter
       */
      virtual jiba::rvec Derivative(const jiba::rvec &FullModel,
          const jiba::rvec &Derivative) const = 0;
      };

    //! This is the simplest transformation, the generalized and physical parameters are identical
    class ModelCopyTransform : public GeneralModelTransform
      {
    public:
      ModelCopyTransform()
        {
        }
      virtual ~ModelCopyTransform()
        {
        }
      virtual jiba::rvec GeneralizedToPhysical(const jiba::rvec &FullModel) const
        {
          return FullModel;
        }
      virtual jiba::rvec PhysicalToGeneralized(const jiba::rvec &FullModel) const
      {
    	return FullModel;
      }
      virtual jiba::rvec Derivative(const jiba::rvec &FullModel,
                const jiba::rvec &Derivative) const
        {
          return Derivative;
        }
      };
    //! A class to handle the distribution of model parameters from the "generalized" parameters of the joint inversion to the concrete parameters for each method
    /*! For the joint inversion we want to have the flexibility to invert for different type of model parameters. For example we can invert
     * for slowness and calculate density and conductivity from it like in the old code, invert for rock physics parameters or use the actual
     * physical parameters for each method.
     *
     * This class handles the translation from these arbitrary inversion parameters to the physical parameters
     * that are required for the forward calculation. It stores a number of function objects, one for each objective
     * function in the inversion that perform the concrete transformation.
     *
     */
    class ModelDistributor
      {
    public:
      ModelDistributor():
        Transformers()
        {
        }
      virtual ~ModelDistributor()
        {
        }
    private:
      //This vector stores a shared pointer to a Transformer object for each objective function
      std::vector<boost::shared_ptr<GeneralModelTransform> > Transformers;
    public:
      //! Add a transformation object
      /*! For each objective function we need to add a transformation object that
       * performs the translation between generalized and physical model parameters.
       * @param Trans A shared pointer to a transformation object
       */
      void AddTransformer(boost::shared_ptr<GeneralModelTransform> Trans)
        {
          Transformers.push_back(Trans);
        }
      //! Perform a transformation using the object stored at index
      /*! In the joint inversion we will have to make sure that we add
       * the transformation objects in the same order as the objective functions.
       * Then this we can translate the model parameters by simply calling this function
       * with the index of the objective function.
       * @param FullModel The generalized model with all parameters
       * @param TransIndex The index of the transformation object
       * @return The physical parameters using the selected transformation
       */
      jiba::rvec operator()(const jiba::rvec &FullModel, const size_t TransIndex)
        {
          assert(TransIndex < Transformers.size());
          return Transformers.at(TransIndex)->GeneralizedToPhysical(FullModel);
        }
      jiba::rvec TransformGradient(const jiba::rvec &FullModel,
          const jiba::rvec &RawGradient, const size_t TransIndex)
        {
          assert(TransIndex < Transformers.size());
          return Transformers.at(TransIndex)->Derivative(FullModel, RawGradient);
        }
      };

  /* @} */
  }

#endif /* MODELDISTRIBUTOR_H_ */
