//============================================================================
// Name        : ResTrans.h
// Author      : Jun 2, 2010
// Version     :
// Copyright   : 2010, mmoorkamp
//============================================================================

#include "../Global/Jif3DGlobal.h"
#include "../Inversion/ModelTransforms.h"

#ifndef RESTRANS_H_
#define RESTRANS_H_
namespace jif3D
  {
    //! The transformation used for calculating velocity and density from resistivity for the synthetic tests developed in Durham
    class J3DEXPORT Durham1DTrans: public jif3D::GeneralModelTransform
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
      //! An object transforming generalized parameters to conductivity before we transform conductivity to velocity
      boost::shared_ptr<GeneralModelTransform> ResTransform;
    public:
      //! Transform the normalized model parameters back to physical parameters
      virtual jif3D::rvec GeneralizedToPhysical(const jif3D::rvec &FullModel) const override
        {
          jif3D::rvec Cond(ResTransform->GeneralizedToPhysical(FullModel));
          jif3D::rvec Output(FullModel.size());
          for (size_t i = 0; i < FullModel.size(); ++i)
            {
              double lres = std::log10(1.0 / Cond(i));
              Output(i) = scale * std::pow(10.0, a + b * lres + c * lres * lres);
            }
          return Output;
        }
      //! Transform the physical model parameters to generalized model parameters
      virtual jif3D::rvec PhysicalToGeneralized(const jif3D::rvec &FullModel) const override
        {

          jif3D::rvec Output(FullModel.size());
          for (size_t i = 0; i < FullModel.size(); ++i)
            {
              double res = (-b + sqrt(b * b + 4.0 * c * (std::log10(FullModel(i)) - a)))
                  / (2.0 * c);
              Output(i) = std::pow(10.0, -res / scale);
            }
          return ResTransform->PhysicalToGeneralized(Output);
        }
      //! Transform the derivative with respect to the physical parameters to normalized parameters
      virtual jif3D::rvec Derivative(const jif3D::rvec &FullModel,
          const jif3D::rvec &Derivative) const override
        {
          jif3D::rvec Res(ResTransform->GeneralizedToPhysical(FullModel));
          jif3D::rvec ResDeriv(ResTransform->Derivative(FullModel, Derivative));
          jif3D::rvec Output(FullModel.size());
          for (size_t i = 0; i < FullModel.size(); ++i)
            {
              Output(i) = std::pow(10.0,
                  a - b * std::log10(Res(i)) + c * std::pow(std::log10(Res(i)), 2))
                  * (2.0 * c * std::log10(Res(i)) / Res(i) - b / Res(i)) * ResDeriv(i)
                  * scale;
            }
          return Output;
        }
      Durham1DTrans(boost::shared_ptr<GeneralModelTransform> ResTrans, const double f1,
          const double f2, const double f3, const double f4) :
          a(f1), b(f2), c(f3), scale(f4), ResTransform(ResTrans)
        {

        }
      };
  }
#endif /* RESTRANS_H_ */
