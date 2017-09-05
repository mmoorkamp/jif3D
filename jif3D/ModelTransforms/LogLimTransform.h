/*
 * LogLimTransform.h
 *
 *  Created on: 13.02.2013
 *      Author:  mm489
 */

#ifndef LOGLIMTRANSFORM_H_
#define LOGLIMTRANSFORM_H_

#include "../Global/Serialization.h"
#include "../Global/Jif3DGlobal.h"
#include "GeneralModelTransform.h"
#include <boost/math/special_functions/atanh.hpp>

namespace jif3D
  {
    /** \addtogroup inversion General routines for inversion */
    /* @{ */

    //! Change an unconstrained optimization problem to a constrained optimization problem through a logarithmic transformation with limits
    /*! We often want to constrain the range of possible values in our optimization problem between an upper and
     * a lower limit. Instead of using a constrained optimization method, we use this transformation
     * \f$ m^{\star} =  \log(m - m_{min}) - \log(m_{max} - m) \f$ between the generalized
     * model parameters \f$ m \f$ and the generalized model parameters \f$ m^{\star} \f$. The generalized model
     * parameters can then vary over the whole range, while the physical model parameters will always been between
     * the two bounds. Note that with this transform the maximum and minimum values will be mapped to infinity
     * so in practice \f$ m_min < m < m_max \f$.
     */
    class J3DEXPORT LogLimTransform: public jif3D::GeneralModelTransform
      {
    private:
      //! The minimum value for the physical model parameters
      double min;
      //! The maximum value for the physical model parameters
      double max;
      friend class access;
      //! Provide serialization to be able to store objects and, more importantly for simpler MPI parallelization
      template<class Archive>
      void serialize(Archive & ar, const unsigned int version)
        {
          ar & base_object<GeneralModelTransform>(*this);
          ar & min;
          ar & max;
        }
    public:
      //! We setup a clone function to have a virtual constructor and create polymorphic copies
      virtual LogLimTransform* clone() const override
        {
          return new LogLimTransform(*this);
        }
      //! Transform the normalized model parameters back to physical parameters
      virtual jif3D::rvec GeneralizedToPhysical(const jif3D::rvec &FullModel) const override
        {
          jif3D::rvec Output(FullModel.size());
          std::transform(FullModel.begin(), FullModel.end(), Output.begin(),
              [this] (double fm)
                {
                  return (max * std::exp(fm) + min)/(std::exp(fm) + 1.0);
                });
          return Output;
        }
      //! Transform the physical model parameters to generalized model parameters
      virtual jif3D::rvec PhysicalToGeneralized(const jif3D::rvec &FullModel) const override
        {
          jif3D::rvec Output(FullModel.size());
          std::transform(FullModel.begin(), FullModel.end(), Output.begin(),
              [this](double fm)
                {
                  return std::log(fm - min) - std::log(max - fm);
                });

          return Output;
        }
      //! Transform the derivative with respect to the physical parameters to normalized parameters
      virtual jif3D::rvec Derivative(const jif3D::rvec &FullModel,
          const jif3D::rvec &Derivative) const override
        {

          jif3D::rvec Output(FullModel.size());
          std::transform(FullModel.begin(), FullModel.end(), Derivative.begin(),
              Output.begin(), [this](double fm, double der)
                { return (max - min) * std::exp(fm)/ std::pow(1.0 + std::exp(fm), 2)
                  * der;});
          return Output;
        }
      //! The constructor need the minimum and maximum physical model parameter
      /*! When we construct the transformation object we need to specify the
       * upper and lower bound for the physical model parameters.
       * @param minval The minimum value for the physical model parameters.
       * @param maxval The maximum value for the physical model parameters.
       */
      LogLimTransform(const double minval = 1.0, const double maxval = 5.0) :
          min(minval), max(maxval)
        {
        }
      virtual ~LogLimTransform()
        {
        }
      };
  /* @} */
  }

#endif /* LOGLIMTRANSFORM_H_ */
