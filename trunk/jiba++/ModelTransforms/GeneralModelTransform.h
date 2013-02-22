/*
 * GeneralModelTransform.h
 *
 *  Created on: 13.02.2013
 *      Author:  mm489
 */

#ifndef GENERALMODELTRANSFORM_H_
#define GENERALMODELTRANSFORM_H_

#include <boost/serialization/serialization.hpp>
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
    private:
      friend class boost::serialization::access;
      //! Provide serialization to be able to store objects and, more importantly for simpler MPI parallelization
      template<class Archive>
      void serialize(Archive & ar, const unsigned int version)
        {

        }
    public:
      //! We setup a clone function to have a virtual constructor and create polymorphic copies
      virtual GeneralModelTransform* clone() const = 0;
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
  /* @} */
  }

#endif /* GENERALMODELTRANSFORM_H_ */
