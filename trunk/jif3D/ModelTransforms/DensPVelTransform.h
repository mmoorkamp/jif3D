/*
 * DensPVelTransform.h
 *
 *  Created on: 28.04.2020
 *      Author:  mm489
 */

#ifndef DENSPVELTRANSFORM_H_
#define DENSPVELTRANSFORM_H_

#include "../Global/Serialization.h"
#include "../Global/Jif3DGlobal.h"
#include "GeneralModelTransform.h"
#include "../ModelBase/ThreeDModelBase.h"
#include <boost/shared_ptr.hpp>

namespace jif3D
  {
    /** \addtogroup inversion General routines for inversion */
    /* @{ */

    //! Transform S-velocity and density anomaly values to S-Velocity, Absolute Density and P-Velocity
    /*!
     */
    class J3DEXPORT DensPVelTransform: public jif3D::GeneralModelTransform
      {
    private:
      //! the vp/vs ratio used to calculate P-wave velocities from vs
      double vpvs;

      friend class access;
      //! Provide serialization to be able to store objects and, more importantly for simpler MPI parallelization
      template<class Archive>
      void serialize(Archive & ar, const unsigned int version)
        {
          ar & base_object<GeneralModelTransform>(*this);
          ar & vpvs;
        }

    public:
      //! We setup a clone function to have a virtual constructor and create polymorphic copies
      virtual DensPVelTransform* clone() const override
        {
          return new DensPVelTransform(*this);
        }
      //! Transform S-velocity and density anomaly values to S-Velocity, P-Velocity and  density anomaly
      virtual jif3D::rvec GeneralizedToPhysical(const jif3D::rvec &FullModel) const
      override
        {
          const size_t nparms = FullModel.size() / 2;
          jif3D::rvec Output(nparms * 3, 0.0);
          std::copy(FullModel.begin(), FullModel.begin() + nparms, Output.begin());
          std::transform(FullModel.begin(), FullModel.begin() + nparms,
              Output.begin() + nparms, [this](double val)
                { return val*vpvs;});
          std::copy(FullModel.begin() + nparms, FullModel.end(),
              Output.begin() + 2 * nparms);

          return Output;
        }
      //! Transform from S-velocity, P-velocity and density to S-velocity and relative density
      virtual jif3D::rvec PhysicalToGeneralized(const jif3D::rvec &FullModel) const
      override
        {
          const size_t nparms = FullModel.size() / 3;
          jif3D::rvec Output(nparms * 2, 0.0);
          std::copy(FullModel.begin(), FullModel.begin() + nparms, Output.begin());
          std::copy(FullModel.begin() + 2 * nparms, FullModel.end(),
              Output.begin() + nparms);

          return Output;
        }
      //! Transform the derivative with respect to the S-velocity, P-velocity and Density to S-vel and density, we ignore vp as it is calculated from vs
      virtual jif3D::rvec Derivative(const jif3D::rvec &FullModel,
          const jif3D::rvec &Derivative) const override
        {
          const size_t nparms = Derivative.size() / 3;
          jif3D::rvec Output(nparms * 2, 0.0);
          std::copy(Derivative.begin(), Derivative.begin() + nparms, Output.begin());
          std::copy(Derivative.begin() + 2* nparms, Derivative.end(),
              Output.begin() + nparms);
          return Output;

        }
      //! The constructor needs a pointer to an object that gives slowness
      /*!
       */
      DensPVelTransform(const std::vector<double> &bg_dens, double vpvs_ = std::sqrt(3.0)) :
          vpvs(vpvs_)
        {
        }
      virtual ~DensPVelTransform()
        {
        }
      };
  /* @} */
  }

#endif /* DENSPVELTRANSFORM_H_ */
