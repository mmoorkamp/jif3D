/*
 * ModelCopyTransform.h
 *
 *  Created on: 13.02.2013
 *      Author:  mm489
 */

#ifndef MODELCOPYTRANSFORM_H_
#define MODELCOPYTRANSFORM_H_

#include <boost/serialization/serialization.hpp>
#include <boost/serialization/base_object.hpp>
#include "../Global/VecMat.h"

namespace jiba
  {
    /** \addtogroup inversion General routines for inversion */
    /* @{ */

    //! This is the simplest transformation, the generalized and physical parameters are identical
    class ModelCopyTransform: public GeneralModelTransform
      {
    private:
      friend class boost::serialization::access;
      //! Provide serialization to be able to store objects and, more importantly for simpler MPI parallelization
      template<class Archive>
      void serialize(Archive & ar, const unsigned int version)
        {
          ar & boost::serialization::base_object<GeneralModelTransform>(*this);
        }
    public:
      //! We setup a clone function to have a virtual constructor and create polymorphic copies
      virtual ModelCopyTransform* clone() const
        {
          return new ModelCopyTransform(*this);
        }
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
    /* @} */
  }

#endif /* MODELCOPYTRANSFORM_H_ */
