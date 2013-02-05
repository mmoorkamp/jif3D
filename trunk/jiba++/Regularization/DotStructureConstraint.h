/*
 * DotStructureConstraint.h
 *
 *  Created on: Feb 5, 2013
 *      Author: mmoorkamp
 */

#ifndef DOTSTRUCTURECONSTRAINT_H_
#define DOTSTRUCTURECONSTRAINT_H_

#include <boost/serialization/serialization.hpp>
#include <boost/serialization/base_object.hpp>
#include "ObjectiveFunction.h"
#include "GradientRegularization.h"

namespace jiba
{

	class DotStructureConstraint: public jiba::ObjectiveFunction
	{
	private:
		//! The object to calculate the spatial gradient for the first model
		GradientRegularization FirstGradient;
		//! The object to calculate the spatial gradient for the second model
		GradientRegularization SecondGradient;
		//! The model geometry for the two models, this is always identical to the geometry in GradientRegularization objects
		ThreeDModelBase ModelGeometry;
		friend class boost::serialization::access;
		//! Provide serialization to be able to store objects and, more importantly for simpler MPI parallelization
		template<class Archive>
		void serialize(Archive & ar, const unsigned int version)
		{
			ar & boost::serialization::base_object<ObjectiveFunction>(*this);

		}
	public:
		//! The implementation of the objective function calculation
		virtual void
		ImplDataDifference(const jiba::rvec &Model, jiba::rvec &Diff);
		//! The gradient of the objective function with respect to the model parameters
		virtual jiba::rvec ImplGradient(const jiba::rvec &Model,
				const jiba::rvec &Diff);
		explicit DotStructureConstraint(const jiba::ThreeDModelBase &Geometry) :
				FirstGradient(Geometry, 0.0), SecondGradient(Geometry, 0.0), ModelGeometry(
						Geometry)
		{

		}
		virtual ~DotStructureConstraint();
	};

} /* namespace jiba */
#endif /* DOTSTRUCTURECONSTRAINT_H_ */
