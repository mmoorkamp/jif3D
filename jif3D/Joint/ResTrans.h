//============================================================================
// Name        : ResTrans.h
// Author      : Jun 2, 2010
// Version     : 
// Copyright   : 2010, mmoorkamp
//============================================================================

#include "../Inversion/ModelTransforms.h"

#ifndef RESTRANS_H_
#define RESTRANS_H_
namespace jiba
  {
    //! The transformation used for calculating velocity and density from resistivity for the synthetic tests developed in Durham
    class Durham1DTrans: public jiba::GeneralModelTransform
      {
    private:
      //! The constant in the conversion equation
      const double a;
      //! The factor for the linear part in the conversion equation
      const double b;
      //! The factor for the square part in the conversion equation
      const double c;
      //! The final scale factor to scale the output by
      const double scale;
      boost::shared_ptr<GeneralModelTransform> ResTransform;
    public:
      //! Transform the normalized model parameters back to physical parameters
      virtual jiba::rvec GeneralizedToPhysical(const jiba::rvec &FullModel) const
        {
          jiba::rvec Cond(ResTransform->GeneralizedToPhysical(FullModel));
          jiba::rvec Output(FullModel.size());
          for (size_t i = 0; i < FullModel.size(); ++i)
            {
              double lres = std::log10(1.0 / Cond(i));
              Output(i) = scale
                  * std::pow(10.0, a + b * lres + c * lres * lres);
            }
          return Output;
        }
      //! Transform the physical model parameters to generalized model parameters
      virtual jiba::rvec PhysicalToGeneralized(const jiba::rvec &FullModel) const
        {

          jiba::rvec Output(FullModel.size());
          for (size_t i = 0; i < FullModel.size(); ++i)
            {
              double res = (-b + sqrt(b * b + 4.0 * c * (std::log10(
                  FullModel(i)) - a))) / (2.0 * c);
              Output(i) = std::pow(10.0, -res / scale);
            }
          return ResTransform->PhysicalToGeneralized(Output);
        }
      //! Transform the derivative with respect to the physical parameters to normalized parameters
      virtual jiba::rvec Derivative(const jiba::rvec &FullModel,
          const jiba::rvec &Derivative) const
        {
          jiba::rvec Res(ResTransform->GeneralizedToPhysical(FullModel));
          jiba::rvec ResDeriv(ResTransform->Derivative(FullModel, Derivative));
          jiba::rvec Output(FullModel.size());
          for (size_t i = 0; i < FullModel.size(); ++i)
            {
              Output(i) = std::pow(10.0, a - b * std::log10(Res(i)) + c
                  * std::pow(std::log10(Res(i)), 2)) * (2.0 * c * std::log10(
                  Res(i)) / Res(i) - b / Res(i)) * ResDeriv(i) * scale;
            }
          return Output;
        }
      Durham1DTrans(boost::shared_ptr<GeneralModelTransform> ResTrans,
          const double f1, const double f2, const double f3, const double f4) :
        a(f1), b(f2), c(f3), scale(f4), ResTransform(ResTrans)
        {

        }
      };
  }
#endif /* RESTRANS_H_ */
