/*
 * DotStructureConstraint.h
 *
 *  Created on: Feb 5, 2013
 *      Author: mmoorkamp
 */

#ifndef DOTSTRUCTURECONSTRAINT_H_
#define DOTSTRUCTURECONSTRAINT_H_

#include "../Global/Jif3DGlobal.h"
#include "../Global/Serialization.h"
#include "../Inversion/ObjectiveFunction.h"
#include "GradientRegularization.h"

namespace jif3D
{

	class J3DEXPORT DotStructureConstraint: public jif3D::ObjectiveFunction
	{
	private:
		//! The object to calculate the spatial gradient for the first model
		GradientRegularization FirstGradient;
		//! The object to calculate the spatial gradient for the second model
		GradientRegularization SecondGradient;
		//! The model geometry for the two models, this is always identical to the geometry in GradientRegularization objects
		ThreeDModelBase ModelGeometry;
		friend class access;
		//! Provide serialization to be able to store objects and, more importantly for simpler MPI parallelization
		template<class Archive>
		void serialize(Archive & ar, const unsigned int version)
		{
			ar & base_object<ObjectiveFunction>(*this);

		}
	public:
    //! The clone function provides a virtual constructor
    virtual DotStructureConstraint *clone() const
      {
        return new DotStructureConstraint(*this);
      }
		//! The implementation of the objective function calculation
		virtual void
		ImplDataDifference(const jif3D::rvec &Model, jif3D::rvec &Diff);
		//! The gradient of the objective function with respect to the model parameters
		virtual jif3D::rvec ImplGradient(const jif3D::rvec &Model,
				const jif3D::rvec &Diff);
		explicit DotStructureConstraint(const jif3D::ThreeDModelBase &Geometry) :
				FirstGradient(Geometry, 0.0), SecondGradient(Geometry, 0.0), ModelGeometry(
						Geometry)
		{

		}
		virtual ~DotStructureConstraint();
	};

} /* namespace jif3D */
#endif /* DOTSTRUCTURECONSTRAINT_H_ */
