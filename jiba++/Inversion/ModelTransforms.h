//============================================================================
// Name        : ModelTransforms.h
// Author      : Apr 15, 2009
// Version     :
// Copyright   : 2009, mmoorkamp
//============================================================================

#ifndef MODELTRANSFORMS_H_
#define MODELTRANSFORMS_H_

#include "ModelDistributor.h"

/*! \file ModelTransforms.h
 * This file contains various classes to transform model parameters within an inversion.
 * These are used to either make the inversion more well behaved or to calculate
 * one physical quantity from another, like density from slowness in LogDensityTransform.
 */

namespace jiba
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
    class NormalizeTransform: public jiba::GeneralModelTransform
      {
    private:
      // the Reference model
      const jiba::rvec Reference;
    public:
      //! Transform the normalized model parameters back to physical parameters
      virtual jiba::rvec GeneralizedToPhysical(const jiba::rvec &FullModel) const
        {
          assert(FullModel.size() == Reference.size());
          return ublas::element_prod(FullModel, Reference);
        }
      virtual jiba::rvec PhysicalToGeneralized(const jiba::rvec &FullModel) const
        {
          assert(FullModel.size() == Reference.size());
          return ublas::element_div(FullModel, Reference);
        }
      //! Transform the derivative with respect to the physical parameters to normalized parameters
      virtual jiba::rvec Derivative(const jiba::rvec &FullModel,
          const jiba::rvec &Derivative) const
        {
          assert(FullModel.size() == Reference.size());
          assert(FullModel.size() == Derivative.size());
          return ublas::element_prod(Derivative, Reference);
        }
      //! The constructor needs the reference model, this has to have the same size as the inversion model
      NormalizeTransform(const jiba::rvec &Ref) :
        Reference(Ref)
        {
        }
      virtual ~NormalizeTransform()
        {
        }
      };
    //! Transform normalized logarithmic parameters
    /*! This transform is used for model parameters of the form \f$m^{\star} = \ln m/m_0\f$.
     * Here \f$m_0\f$ is the reference mdoel, e.g. the starting model of the inversion.
     * The advantage over simple normalization is that we enforce the physical parameters
     * to be positive and that we can capture large variations in physical parameters
     * in a relatively small range of inversion parameters.
     */
    class LogTransform: public jiba::GeneralModelTransform
      {
    private:
      const jiba::rvec Reference;
    public:
      //! Transform the normalized model parameters back to physical parameters
      virtual jiba::rvec GeneralizedToPhysical(const jiba::rvec &FullModel) const
        {
          assert(FullModel.size() == Reference.size());
          jiba::rvec Output(FullModel.size());
          for (size_t i = 0; i < FullModel.size(); ++i)
            Output( i) = exp(FullModel(i)) * Reference(i);
          return Output;
        }
      virtual jiba::rvec PhysicalToGeneralized(const jiba::rvec &FullModel) const
        {
          assert(FullModel.size() == Reference.size());
          jiba::rvec Output(FullModel.size());
          for (size_t i = 0; i < FullModel.size(); ++i)
            Output( i) = log(FullModel(i) / Reference(i));
          return Output;
        }
      //! Transform the derivative with respect to the physical parameters to normalized parameters
      virtual jiba::rvec Derivative(const jiba::rvec &FullModel,
          const jiba::rvec &Derivative) const
        {

          jiba::rvec Output(FullModel.size());
          for (size_t i = 0; i < FullModel.size(); ++i)
            {
              Output( i) = Reference(i) * exp(FullModel(i)) * Derivative(i);
            }
          return Output;
        }
      LogTransform(const jiba::rvec &Ref) :
        Reference(Ref)
        {
        }
      virtual ~LogTransform()
        {
        }
      };

    class TanhTransform: public jiba::GeneralModelTransform
          {
        private:
        	const double min;
        	const double max;
          const jiba::rvec Reference;
        public:
          //! Transform the normalized model parameters back to physical parameters
          virtual jiba::rvec GeneralizedToPhysical(const jiba::rvec &FullModel) const
            {
              assert(FullModel.size() == Reference.size());
              jiba::rvec Output(FullModel.size());
              for (size_t i = 0; i < FullModel.size(); ++i)
                Output( i) = min + (1.0 + tanh(FullModel(i)))/2.0 * (max-min);
              return Output;
            }
          virtual jiba::rvec PhysicalToGeneralized(const jiba::rvec &FullModel) const
            {
              assert(FullModel.size() == Reference.size());
              jiba::rvec Output(FullModel.size());
              for (size_t i = 0; i < FullModel.size(); ++i)
              {
            	  const double argument = 2.0 * (FullModel(i)- min)/(max-min)  - 1;
                Output( i) = atanh(argument);
              }
              return Output;
            }
          //! Transform the derivative with respect to the physical parameters to normalized parameters
          virtual jiba::rvec Derivative(const jiba::rvec &FullModel,
              const jiba::rvec &Derivative) const
            {

              jiba::rvec Output(FullModel.size());
              for (size_t i = 0; i < FullModel.size(); ++i)
                {
                  Output( i) = (max-min)/(2.0 * pow(cosh(FullModel(i)),2)) * Derivative(i);
                }
              return Output;
            }
          TanhTransform(const jiba::rvec &Ref, const double minval = 1.0, const double maxval = 5.0) :
            min(minval),max(maxval),Reference(Ref)
            {
            }
          virtual ~TanhTransform()
            {
            }
          };

    class TanhDensityTransform: public jiba::TanhTransform
{
public:
	virtual jiba::rvec GeneralizedToPhysical(const jiba::rvec &FullModel) const
	{
		jiba::rvec Slowness(TanhTransform::GeneralizedToPhysical(FullModel));
		jiba::rvec Output(FullModel.size());
		for (size_t i = 0; i < FullModel.size(); ++i)
		{
			Output(i) = (1.0 / Slowness(i) + 8500.0) / 5000.0;
		}
		return Output;
	}
	virtual jiba::rvec PhysicalToGeneralized(const jiba::rvec &FullModel) const
	{

		jiba::rvec Output(FullModel.size());
				for (size_t i = 0; i < FullModel.size(); ++i)
				{
					Output(i) = 1.0/(5000.0 * FullModel(i) - 8500.0);
				}
				return TanhTransform::PhysicalToGeneralized(Output);
	}
	virtual jiba::rvec Derivative(const jiba::rvec &FullModel,
			const jiba::rvec &Derivative) const
	{
		jiba::rvec Slowness(TanhTransform::GeneralizedToPhysical(FullModel));
		jiba::rvec SlowDeriv(TanhTransform::Derivative(FullModel, Derivative));
		jiba::rvec Output(FullModel.size());
			for (size_t i = 0; i < FullModel.size(); ++i)
			{
				Output(i) = -1.0 / (Slowness(i) * Slowness(i))  / 5000.0 * SlowDeriv(i);
			}
			return Output;
	}
	TanhDensityTransform(const jiba::rvec &Ref, const double minval = 1.0,
			const double maxval = 5.0) :
		TanhTransform(Ref, minval, maxval)
	{
	}
	virtual ~TanhDensityTransform()
	{
	}
    };


    class LogDensityTransform: public jiba::GeneralModelTransform
      {
    private:
      const jiba::rvec Reference;
    public:
      virtual jiba::rvec GeneralizedToPhysical(const jiba::rvec &FullModel) const
        {
          jiba::rvec Output(FullModel.size());
          for (size_t i = 0; i < FullModel.size(); ++i)
            Output( i) = (1. / (exp(FullModel(i)) * Reference(i)) + 8500.0)
                / 5000.0;
          return Output;
        }
      virtual jiba::rvec PhysicalToGeneralized(const jiba::rvec &FullModel) const
        {
          jiba::rvec Output(FullModel.size());
          for (size_t i = 0; i < FullModel.size(); ++i)
            Output( i) = log(1.0 / (Reference(i) * (5000.0 * FullModel(i)
                - 8500.0)));
          return Output;
        }
      virtual jiba::rvec Derivative(const jiba::rvec &FullModel,
          const jiba::rvec &Derivative) const
        {
          const double factor = -1.0 / 5000.0;
          jiba::rvec Output(FullModel.size());
          for (size_t i = 0; i < FullModel.size(); ++i)
            {
              Output( i) = factor / (Reference(i) * exp(FullModel(i)))
                  * Derivative(i);
            }
          return Output;
        }
      LogDensityTransform(const jiba::rvec &Ref) :
        Reference(Ref)
        {
        }
      virtual ~LogDensityTransform()
        {
        }
      };

    class VelTransform: public jiba::GeneralModelTransform
      {
    private:
      const jiba::rvec Reference;
    public:
      virtual jiba::rvec GeneralizedToPhysical(const jiba::rvec &FullModel) const
        {
          jiba::rvec Output(FullModel.size());
          for (size_t i = 0; i < FullModel.size(); ++i)
            Output( i) = 1.0 / (FullModel(i) * Reference(i));
          return Output;
        }
      virtual jiba::rvec PhysicalToGeneralized(const jiba::rvec &FullModel) const
        {
          jiba::rvec Output(FullModel.size());
          for (size_t i = 0; i < FullModel.size(); ++i)
            Output( i) = 1.0 / (FullModel(i) * Reference(i));
          return Output;
        }
      virtual jiba::rvec Derivative(const jiba::rvec &FullModel,
          const jiba::rvec &Derivative) const
        {
          jiba::rvec Output(FullModel.size());
          for (size_t i = 0; i < FullModel.size(); ++i)
            Output( i) = -1.0 / (Reference(i) * FullModel(i) * FullModel(i))
                * Derivative(i);
          return Output;
        }
      VelTransform(const jiba::rvec &Ref) :
        Reference(Ref)
        {
        }
      virtual ~VelTransform()
        {
        }
      };

    class VelDensTransform: public jiba::GeneralModelTransform
      {
    private:
      jiba::rvec Reference;
    public:
      virtual jiba::rvec GeneralizedToPhysical(const jiba::rvec &FullModel) const
        {
          jiba::rvec Output(FullModel.size());
          for (size_t i = 0; i < FullModel.size(); ++i)
            {
              Output( i) = (Reference(i) * FullModel(i) + 8500.0) / 5000.0;
            }
          return Output;
        }
      virtual jiba::rvec PhysicalToGeneralized(const jiba::rvec &FullModel) const
        {
          jiba::rvec Output(FullModel.size());
          for (size_t i = 0; i < FullModel.size(); ++i)
            {
              Output( i) = (5000.0 * FullModel(i) - 8500) / Reference(i);
            }
          return Output;
        }
      virtual jiba::rvec Derivative(const jiba::rvec &FullModel,
          const jiba::rvec &Derivative) const
        {
          jiba::rvec Output(FullModel.size());
          for (size_t i = 0; i < FullModel.size(); ++i)
            {
              Output( i) = Reference(i) / 5000.0 * Derivative(i);
            }
          return Output;
        }
      VelDensTransform(const jiba::rvec &Ref) :
        Reference(Ref)
        {
        }
      virtual ~VelDensTransform()
        {
        }
      };
  /* @} */
  }
#endif /* MODELTRANSFORMS_H_ */
