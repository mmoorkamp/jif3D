//============================================================================
// Name        : ResTrans.h
// Author      : Jun 2, 2010
// Version     : 
// Copyright   : 2010, mmoorkamp
//============================================================================

#include "../Inversion/ModelTransforms.h"
#include <boost/math/special_functions/pow.hpp>
#include <cmath>

#ifndef RESTRANS_H_
#define RESTRANS_H_
namespace jiba
  {
    //! The transformation used for calculating velocity and density from resistivity for the synthetic tests developed in Durham
    class SlowCondTrans: public jiba::GeneralModelTransform
      {
    private:
      //! The constant in the conversion equation
      const double a;
      //! The factor for the linear part in the conversion equation
      const double b;
      //! The factor for the square part in the conversion equation
      const double c;

      boost::shared_ptr<GeneralModelTransform> SlowTransform;
    public:
      //! We setup a clone function to have a virtual constructor and create polymorphic copies
      virtual SlowCondTrans* clone() const
        {
          return new SlowCondTrans(*this);
        }
      //! Transform the normalized model parameters back to physical parameters
      virtual jiba::rvec GeneralizedToPhysical(const jiba::rvec &FullModel) const
        {
          using boost::math::pow;
          jiba::rvec Slow(SlowTransform->GeneralizedToPhysical(FullModel));
          jiba::rvec Output(FullModel.size());
          for (size_t i = 0; i < FullModel.size(); ++i)
            {
              //std::cout << 1.0 / Slow(i) << " ";
              Output(i) = std::exp(-a - b / Slow(i) - c / pow<2>(Slow(i)));

            }
          return Output;
        }
      //! Transform the physical model parameters to generalized model parameters
      virtual jiba::rvec PhysicalToGeneralized(const jiba::rvec &FullModel) const
        {

          jiba::rvec Output(FullModel.size());
          for (size_t i = 0; i < FullModel.size(); ++i)
            {
              Output(i) = (2.0 * c)
                  / (-b + sqrt(b * b - 4.0 * c * (std::log(FullModel(i)) + a)));
            }
          return SlowTransform->PhysicalToGeneralized(Output);
        }
      //! Transform the derivative with respect to the physical parameters to normalized parameters
      virtual jiba::rvec Derivative(const jiba::rvec &FullModel,
          const jiba::rvec &Derivative) const
        {
          jiba::rvec Slow(SlowTransform->GeneralizedToPhysical(FullModel));
          jiba::rvec SlowDeriv(SlowTransform->Derivative(FullModel, Derivative));
          jiba::rvec Output(FullModel.size());
          for (size_t i = 0; i < FullModel.size(); ++i)
            {
              using boost::math::pow;
              Output(i) = std::exp(-a - b / Slow(i) - c / pow<2>(Slow(i)))
                  * (2.0 * c / pow<3>(Slow(i)) + b / Slow(i)) * SlowDeriv(i);
            }
          return Output;
        }
      SlowCondTrans(boost::shared_ptr<GeneralModelTransform> SlowTrans, const double f1 =
          5.6217, const double f2 = -0.0029885, const double f3 = 4.8008e-7) :
          a(f1), b(f2), c(f3), SlowTransform(SlowTrans)
        {

        }
      };
  }
#endif /* RESTRANS_H_ */
