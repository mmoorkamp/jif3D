/*
 * LogTransform.h
 *
 *  Created on: 13.02.2013
 *      Author:  mm489
 */

#ifndef LOGTRANSFORM_H_
#define LOGTRANSFORM_H_

#include "../Global/Serialization.h"
#include "GeneralModelTransform.h"

namespace jif3D
  {
    /** \addtogroup inversion General routines for inversion */
    /* @{ */

    //! Transform normalized logarithmic parameters
    /*! This transform is used for model parameters of the form \f$m^{\star} = \ln m/m_0\f$.
     * Here \f$m_0\f$ is the reference mdoel, e.g. the starting model of the inversion.
     * The advantage over simple normalization is that we enforce the physical parameters
     * to be positive and that we can capture large variations in physical parameters
     * in a relatively small range of inversion parameters.
     */
    class LogTransform: public jif3D::GeneralModelTransform
      {
    private:
      //! Each model parameter is divided by the reference values before taking the logarithm
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
      virtual LogTransform* clone() const
        {
          return new LogTransform(*this);
        }
      //! Transform the normalized model parameters back to physical parameters
      virtual jif3D::rvec GeneralizedToPhysical(const jif3D::rvec &FullModel) const
        {
          assert(FullModel.size() == Reference.size());
          jif3D::rvec Output(FullModel.size());
          for (size_t i = 0; i < FullModel.size(); ++i)
            Output(i) = std::exp(FullModel(i)) * Reference(i);
          return Output;
        }
      //! Transform the physical model parameters to generalized model parameters
      virtual jif3D::rvec PhysicalToGeneralized(const jif3D::rvec &FullModel) const
        {
          assert(FullModel.size() == Reference.size());
          jif3D::rvec Output(FullModel.size());
          for (size_t i = 0; i < FullModel.size(); ++i)
            Output(i) = std::log(FullModel(i) / Reference(i));
          return Output;
        }
      //! Transform the derivative with respect to the physical parameters to normalized parameters
      virtual jif3D::rvec Derivative(const jif3D::rvec &FullModel,
          const jif3D::rvec &Derivative) const
        {

          jif3D::rvec Output(FullModel.size());
          for (size_t i = 0; i < FullModel.size(); ++i)
            {
              Output(i) = Reference(i) * std::exp(FullModel(i)) * Derivative(i);
            }
          return Output;
        }
      //! The constructor takes a reference model for normalization as argument
      /*! When we construct the transformation object we need the values for \f$m_0\f$ as described in
       * the general description for this class.
       * @param Ref The reference model vector \f$m_0\f$
       */
      LogTransform(const jif3D::rvec &Ref) :
          Reference(Ref)
        {
        }
      virtual ~LogTransform()
        {
        }
      };
  /* @} */
  }

#endif /* LOGTRANSFORM_H_ */
