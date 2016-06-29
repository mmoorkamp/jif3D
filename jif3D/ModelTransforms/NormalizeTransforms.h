/*
 * NormalizeTransforms.h
 *
 *  Created on: 13.02.2013
 *      Author:  mm489
 */

#ifndef NORMALIZETRANSFORMS_H_
#define NORMALIZETRANSFORMS_H_

#include "../Global/Serialization.h"
#include "../Global/Jif3DGlobal.h"
#include "GeneralModelTransform.h"

namespace jif3D
  {
    /** \addtogroup inversion General routines for inversion */
    /* @{ */
    //! Normalize the model parameters by dividing by a reference model.
    /*! This class takes a reference model, e.g. the starting model in the
     * inversion and divides each model parameter by the corresponding value
     * in the reference model. The two models therefore have to have the same length.
     * This makes the model parameters dimensionless and, at least in the beginning, on the
     * order of unity therefore helping to avoid problems with greatly varying magnitudes.
     */
    class J3DEXPORT NormalizeTransform: public jif3D::GeneralModelTransform
      {
    private:
      //! The Reference model we devide the model parameters by
      const jif3D::rvec Reference;
      friend class access;
      //! Provide serialization to be able to store objects and, more importantly for simpler MPI parallelization
      template<class Archive>
      void serialize(Archive & ar, const unsigned int version)
        {
          ar & base_object<GeneralModelTransform>(*this);
          ar & Reference;
        }
    public:
      //! We setup a clone function to have a virtual constructor and create polymorphic copies
      virtual NormalizeTransform* clone() const
        {
          return new NormalizeTransform(*this);
        }
      //! Transform the normalized model parameters back to physical parameters
      virtual jif3D::rvec GeneralizedToPhysical(const jif3D::rvec &FullModel) const
        {
          assert(FullModel.size() == Reference.size());
          return ublas::element_prod(FullModel, Reference);
        }
      //! Transform the physical model parameters to generalized model parameters
      virtual jif3D::rvec PhysicalToGeneralized(const jif3D::rvec &FullModel) const
        {
          assert(FullModel.size() == Reference.size());
          return ublas::element_div(FullModel, Reference);
        }
      //! Transform the derivative with respect to the physical parameters to normalized parameters
      virtual jif3D::rvec Derivative(const jif3D::rvec &FullModel,
          const jif3D::rvec &Derivative) const
        {
          return GeneralizedToPhysical(Derivative);
        }
      //! The constructor needs the reference model, this has to have the same size as the inversion model
      NormalizeTransform(const jif3D::rvec &Ref) :
          Reference(Ref)
        {
        }
      virtual ~NormalizeTransform()
        {
        }
      };
  /* @} */
  }

#endif /* NORMALIZETRANSFORMS_H_ */
