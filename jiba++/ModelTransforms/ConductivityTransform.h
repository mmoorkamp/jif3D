/*
 * ConductivityTransform.h
 *
 *  Created on: 13.02.2013
 *      Author: mm489
 */

#ifndef CONDUCTIVITYTRANSFORM_H_
#define CONDUCTIVITYTRANSFORM_H_

#include <boost/serialization/serialization.hpp>
#include <boost/serialization/base_object.hpp>
#include <boost/shared_ptr.hpp>
#include "GeneralModelTransform.h"


namespace jiba
{
	/** \addtogroup inversion General routines for inversion */
	/* @{ */

//! Transform generalized model parameters to conductivity
	/*! Similarly to DensityTransform, this class transforms
	 * generalized model parameters to conductivity taking an
	 * intermediate step through slowness. This is motivated
	 * by the parametrization used in the SINDRI project.
	 * The equation for the transform is
	 * \f$ \sigma = \exp( -a/s^2 - b/s -c), \f$
	 * where s is slowness and the coefficient a,b,c can be specified
	 * in the constructor.
	 */
	class ConductivityTransform: public jiba::GeneralModelTransform
	{
	private:
		//! The coefficient for the quadratic term of the transform
		const double a;
		//! The coefficient for the linear term  of the transform
		const double b;
		//! The constant term of the transform
		const double c;
		boost::shared_ptr<GeneralModelTransform> SlownessTransform;
		//! An object indicating where to apply the parameter relationship (value 1)
		jiba::ThreeDModelBase RelModel;
		//! The value to use for density where the relationship is not valid
		double replacevalue;
		friend class boost::serialization::access;
		//! Provide serialization to be able to store objects and, more importantly for simpler MPI parallelization
		template<class Archive>
		void serialize(Archive & ar, const unsigned int version)
		{
			ar
					& boost::serialization::base_object<GeneralModelTransform>(
							*this);
			ar & a;
			ar & b;
			ar & c;
			ar & SlownessTransform;
			ar & RelModel;
			ar & replacevalue;
		}
	public:
		//! We setup a clone function to have a virtual constructor and create polymorphic copies
		virtual ConductivityTransform* clone() const
		{
			return new ConductivityTransform(*this);
		}
		//! Transform Generalized parameters in terms of slowness to conductivity using a functional relationship
		virtual jiba::rvec GeneralizedToPhysical(
				const jiba::rvec &FullModel) const
		{
			jiba::rvec Slowness(
					SlownessTransform->GeneralizedToPhysical(FullModel));
			jiba::rvec Conductivity(FullModel.size());
			for (size_t i = 0; i < FullModel.size(); ++i)
			{
				if (RelModel.GetData().data()[i])
				{
				Conductivity(i) =
						std::exp(
								-(a / (Slowness(i) * Slowness(i))
										+ b / Slowness(i) + c));
				}
				else
				{
					Conductivity(i) = replacevalue;
				}
			}
			return Conductivity;
		}
		//! Transform Conductivity to Slowness and then Generalized Parameters
		virtual jiba::rvec PhysicalToGeneralized(
				const jiba::rvec &FullModel) const
		{
			size_t nvals = FullModel.size();
			jiba::rvec Slowness(nvals);
			for (size_t i = 0; i < nvals; ++i)
			{
				if (RelModel.GetData().data()[i])
				{
				double vel = (-b
						+ sqrt(b * b - 4.0 * a * (c + std::log(FullModel(i)))))
						/ (2.0 * a);
				Slowness(i) = 1.0 / vel;
				}
				else
				{
					Slowness(i) = replacevalue;
				}
			}
			return SlownessTransform->PhysicalToGeneralized(Slowness);
		}
		//! Transform the derivative with respect to the physical parameters to normalized parameters
		virtual jiba::rvec Derivative(const jiba::rvec &FullModel,
				const jiba::rvec &Derivative) const
		{
			jiba::rvec Slowness(
					SlownessTransform->GeneralizedToPhysical(FullModel));
			jiba::rvec SlowDeriv(
					SlownessTransform->Derivative(FullModel, Derivative));
			jiba::rvec Gradient(FullModel.size());
			for (size_t i = 0; i < FullModel.size(); ++i)
			{
				if (RelModel.GetData().data()[i])
				{
				Gradient(i) =
						std::exp(
								-(a / (Slowness(i) * Slowness(i))
										+ b / Slowness(i) + c))
								* (2.0 * a * std::pow(Slowness(i), -3)
										+ b / (Slowness(i) * Slowness(i)))
								* SlowDeriv(i);
				}
				else
				{
					Gradient = 0.0;
				}
			}
			return Gradient;
		}
		//! The constructor needs a pointer to an object that gives slowness
		/*! To reduce the amount of code in the main program this class
		 * takes a pointer to a transformation object that it uses to
		 * transform the generalized model parameters to slowness before
		 * then transforming to conductivity. In addition we can specify
		 * the coefficients for the transformation.
		 * We assume a functional relationship of the form
		 * \f$ \rho = \exp( -a/v^2 + b/v +c) \f$.
		 * @param SlowTrans A pointer to an object that gives slowness
		 * @param aval The coefficient for the quadratic term
		 * @param bval The coefficient for the linear term
		 * @param cval The offset in the exponential
		 */
		ConductivityTransform(
				boost::shared_ptr<GeneralModelTransform> SlowTrans,
				const jiba::ThreeDModelBase &RModel, double rvalue,
				double aval = 2.31e-7, double bval = -5.79e-4, double cval =
						0.124) :
				a(aval), b(bval), c(cval), SlownessTransform(SlowTrans), replacevalue(rvalue),
				RelModel(RModel)
		{
		}
		virtual ~ConductivityTransform()
		{
		}
	};
	/* @} */
	}

#endif /* CONDUCTIVITYTRANSFORM_H_ */
