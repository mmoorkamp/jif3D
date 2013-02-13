/*
 * LogTransform.h
 *
 *  Created on: 13.02.2013
 *      Author:  mm489
 */

#ifndef LOGTRANSFORM_H_
#define LOGTRANSFORM_H_

#include <boost/serialization/serialization.hpp>
#include <boost/serialization/base_object.hpp>
#include "GeneralModelTransform.h"

namespace jiba
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
	class LogTransform: public jiba::GeneralModelTransform
	{
	private:
		//! Each model parameter is divided by the reference values before taking the logarithm
		const jiba::rvec Reference;
		friend class boost::serialization::access;
		//! Provide serialization to be able to store objects and, more importantly for simpler MPI parallelization
		template<class Archive>
		void serialize(Archive & ar, const unsigned int version)
		{
			ar
					& boost::serialization::base_object<GeneralModelTransform>(
							*this);
			ar & Reference;
		}
	public:
		//! We setup a clone function to have a virtual constructor and create polymorphic copies
		virtual LogTransform* clone() const
		{
			return new LogTransform(*this);
		}
		//! Transform the normalized model parameters back to physical parameters
		virtual jiba::rvec GeneralizedToPhysical(
				const jiba::rvec &FullModel) const
		{
			assert(FullModel.size() == Reference.size());
			jiba::rvec Output(FullModel.size());
			for (size_t i = 0; i < FullModel.size(); ++i)
				Output(i) = std::exp(FullModel(i)) * Reference(i);
			return Output;
		}
		//! Transform the physical model parameters to generalized model parameters
		virtual jiba::rvec PhysicalToGeneralized(
				const jiba::rvec &FullModel) const
		{
			assert(FullModel.size() == Reference.size());
			jiba::rvec Output(FullModel.size());
			for (size_t i = 0; i < FullModel.size(); ++i)
				Output(i) = std::log(FullModel(i) / Reference(i));
			return Output;
		}
		//! Transform the derivative with respect to the physical parameters to normalized parameters
		virtual jiba::rvec Derivative(const jiba::rvec &FullModel,
				const jiba::rvec &Derivative) const
		{

			jiba::rvec Output(FullModel.size());
			for (size_t i = 0; i < FullModel.size(); ++i)
			{
				Output(i) = Reference(i) * std::exp(FullModel(i))
						* Derivative(i);
			}
			return Output;
		}
		//! The constructor takes a reference model for normalization as argument
		/*! When we construct the transformation object we need the values for \f$m_0\f$ as described in
		 * the general description for this class.
		 * @param Ref The reference model vector \f$m_0\f$
		 */
		LogTransform(const jiba::rvec &Ref) :
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
