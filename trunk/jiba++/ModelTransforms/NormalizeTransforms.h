/*
 * NormalizeTransforms.h
 *
 *  Created on: 13.02.2013
 *      Author:  mm489
 */

#ifndef NORMALIZETRANSFORMS_H_
#define NORMALIZETRANSFORMS_H_

#include <boost/serialization/serialization.hpp>
#include <boost/serialization/base_object.hpp>
#include "GeneralModelTransform.h"

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
		//! The Reference model we devide the model parameters by
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
		virtual NormalizeTransform* clone() const
		{
			return new NormalizeTransform(*this);
		}
		//! Transform the normalized model parameters back to physical parameters
		virtual jiba::rvec GeneralizedToPhysical(
				const jiba::rvec &FullModel) const
		{
			assert(FullModel.size() == Reference.size());
			return ublas::element_prod(FullModel, Reference);
		}
		//! Transform the physical model parameters to generalized model parameters
		virtual jiba::rvec PhysicalToGeneralized(
				const jiba::rvec &FullModel) const
		{
			assert(FullModel.size() == Reference.size());
			return ublas::element_div(FullModel, Reference);
		}
		//! Transform the derivative with respect to the physical parameters to normalized parameters
		virtual jiba::rvec Derivative(const jiba::rvec &FullModel,
				const jiba::rvec &Derivative) const
		{
			return GeneralizedToPhysical(Derivative);
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
	/* @} */
}

#endif /* NORMALIZETRANSFORMS_H_ */
