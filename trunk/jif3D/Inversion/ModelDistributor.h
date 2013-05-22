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
#include <boost/shared_ptr.hpp>
#include <boost/serialization/serialization.hpp>
#include <boost/serialization/shared_ptr.hpp>
#include "../Global/VecMat.h"
#include "../ModelTransforms/GeneralModelTransform.h"

namespace jif3D
  {
    /** \addtogroup inversion General routines for inversion */
    /* @{ */

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
      ModelDistributor() :
          Transformers()
        {
        }
      virtual ~ModelDistributor()
        {
        }
    private:
      //! This vector stores a shared pointer to a Transformer object for each objective function
      std::vector<boost::shared_ptr<GeneralModelTransform> > Transformers;
      friend class boost::serialization::access;
      //! Provide serialization to be able to store objects and, more importantly for simpler MPI parallelization
      template<class Archive>
      void serialize(Archive & ar, const unsigned int version)
        {
          ar & Transformers;
        }
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
      jif3D::rvec operator()(const jif3D::rvec &FullModel, const size_t TransIndex)
        {
          assert(TransIndex < Transformers.size());
          return Transformers.at(TransIndex)->GeneralizedToPhysical(FullModel);
        }
      //! Transform the gradient for a given model using the transformation stored at index
      /*! When we calculate the gradient of the objective function for a given method, we
       * get the gradient with respect to the physical parameters. We can use this function to
       * transform the gradient to the inversion parameters.
       * @param FullModel The generalized model with all parameters
       * @param RawGradient The gradient with respect to the physical parameters
       * @param TransIndex The index of the transformation object
       * @return The gradient with respect to the inversion parameters
       */
      jif3D::rvec TransformGradient(const jif3D::rvec &FullModel,
          const jif3D::rvec &RawGradient, const size_t TransIndex)
        {
          assert(TransIndex < Transformers.size());
          return Transformers.at(TransIndex)->Derivative(FullModel, RawGradient);
        }
      };

  /* @} */
  }

#endif /* MODELDISTRIBUTOR_H_ */
