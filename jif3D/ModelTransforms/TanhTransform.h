/*
 * TanhTransform.h
 *
 *  Created on: 13.02.2013
 *      Author:  mm489
 */

#ifndef TANHTRANSFORM_H_
#define TANHTRANSFORM_H_

#include "../Global/Serialization.h"
#include "../Global/Jif3DGlobal.h"
#include "../Global/VecMat.h"
#include "GeneralModelTransform.h"
#include <boost/math/special_functions/atanh.hpp>

namespace jif3D
  {
    /** \addtogroup inversion General routines for inversion */
    /* @{ */

    //! Change an unconstrained optimization problem to a constrained optimization problem through a tanh transformation
    /*! We often want to constrain the range of possible values in our optimization problem between an upper and
     * a lower limit. Instead of using a constrained optimization method, we use this transformation
     * \f$ m^{\star} =  \mbox{atanh} \left(2.0 * \frac{m - m_{min}} {m_{max} - m_{min}} \right) \f$ between the generalized
     * model parameters \f$ m \f$ and the generalized model parameters \f$ m^{\star} \f$. The generalized model
     * parameters can then vary over the whole range, while the physical model parameters will always been between
     * the two bounds. Note that with this transform the maximum and minimum values will be mapped to infinity
     * so in practice \f$ m_min < m < m_max \f$.
     */
    class J3DEXPORT TanhTransform: public jif3D::GeneralModelTransform
      {
    private:
      //! The minimum value for the physical model parameters
      mutable jif3D::rvec min;
      //! The maximum value for the physical model parameters
      mutable jif3D::rvec max;
      friend class access;
      //! Provide serialization to be able to store objects and, more importantly for simpler MPI parallelization
      template<class Archive>
      void serialize(Archive & ar, const unsigned int version)
        {
          ar & base_object<GeneralModelTransform>(*this);
          ar & min;
          ar & max;
        }
      void CheckMinMax(size_t nelem) const
        {
          if (min.size() != nelem)
            {
              if (min.size() != 1)
                {
                  jif3D::FatalException("Problem in minimum/maximum specification",
                  __FILE__, __LINE__);
                }
              else
                {
                  double minv = min(0);
                  double maxv = max(0);
                  min.resize(nelem);
                  max.resize(nelem);
                  std::fill(min.begin(), min.end(), minv);
                  std::fill(max.begin(), max.end(), maxv);
                }
            }
        }
    public:
      //! We setup a clone function to have a virtual constructor and create polymorphic copies
      virtual TanhTransform* clone() const override
        {
          return new TanhTransform(*this);
        }
      //! Transform the normalized model parameters back to physical parameters
      virtual jif3D::rvec GeneralizedToPhysical(const jif3D::rvec &FullModel) const
      override
        {
          const size_t nelem = FullModel.size();
          CheckMinMax(nelem);
          jif3D::rvec Output(nelem);
          for (size_t i = 0; i < nelem; ++i)
            {
              Output(i) = min(i) + (1.0 + std::tanh(FullModel(i))) / 2.0 * (max(i) - min(i));
            }
          return Output;
        }
      //! Transform the physical model parameters to generalized model parameters
      virtual jif3D::rvec PhysicalToGeneralized(const jif3D::rvec &FullModel) const
      override
        {
          const size_t nelem = FullModel.size();
          CheckMinMax(nelem);
          jif3D::rvec Output(nelem);
          for (size_t i = 0; i < nelem; ++i)
            {
              Output(i) = boost::math::atanh(
                  2.0 * (FullModel(i) - min(i)) / (max(i) - min(i)) - 1.0);
            }

          return Output;
        }
      //! Transform the derivative with respect to the physical parameters to normalized parameters
      virtual jif3D::rvec Derivative(const jif3D::rvec &FullModel,
          const jif3D::rvec &Derivative) const override
        {
          const size_t nelem = FullModel.size();
          CheckMinMax(nelem);
          jif3D::rvec Output(nelem);
          for (size_t i = 0; i < nelem; ++i)
            {
              Output(i) = (max(i) - min(i)) / (2.0 * std::pow(std::cosh(FullModel(i)), 2))
                  * Derivative(i);
            }

          return Output;
        }
      //! The constructor need the minimum and maximum physical model parameter
      /*! When we construct the transformation object we need to specify the
       * upper and lower bound for the physical model parameters.
       * @param minval The minimum value for the physical model parameters.
       * @param maxval The maximum value for the physical model parameters.
       */
      TanhTransform(double minval, double maxval) :
          min(1, minval), max(1, maxval)
        {
        }
      TanhTransform(jif3D::rvec minval, jif3D::rvec maxval) :
          min(minval), max(maxval)
        {
        }
      virtual ~TanhTransform()
        {
        }
      };
  /* @} */
  }

#endif /* TANHTRANSFORM_H_ */
