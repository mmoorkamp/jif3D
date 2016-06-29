/*
 * ModelCopyTransform.h
 *
 *  Created on: 13.02.2013
 *      Author:  mm489
 */

#ifndef MODELCOPYTRANSFORM_H_
#define MODELCOPYTRANSFORM_H_

#include "../Global/Serialization.h"
#include "../Global/Jif3DGlobal.h"
#include "../Global/VecMat.h"
#include "GeneralModelTransform.h"

namespace jif3D
  {
    /** \addtogroup inversion General routines for inversion */
    /* @{ */

    //! This is the simplest transformation, the generalized and physical parameters are identical
    class J3DEXPORT ModelCopyTransform: public GeneralModelTransform
      {
    private:
      friend class access;
      //! Provide serialization to be able to store objects and, more importantly for simpler MPI parallelization
      template<class Archive>
      void serialize(Archive & ar, const unsigned int version)
        {
          ar & base_object<GeneralModelTransform>(*this);
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
      virtual jif3D::rvec GeneralizedToPhysical(const jif3D::rvec &FullModel) const
        {
          return FullModel;
        }
      virtual jif3D::rvec PhysicalToGeneralized(const jif3D::rvec &FullModel) const
        {
          return FullModel;
        }
      virtual jif3D::rvec Derivative(const jif3D::rvec &FullModel,
          const jif3D::rvec &Derivative) const
        {
          return Derivative;
        }
      };
  /* @} */
  }

#endif /* MODELCOPYTRANSFORM_H_ */
