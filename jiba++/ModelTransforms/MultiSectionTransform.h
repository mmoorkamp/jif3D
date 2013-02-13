/*
 * MultiSectionTransform.h
 *
 *  Created on: 13.02.2013
 *      Author:  mm489
 */

#ifndef MULTISECTIONTRANSFORM_H_
#define MULTISECTIONTRANSFORM_H_

#include <numeric>
#include <vector>
#include <boost/serialization/serialization.hpp>
#include <boost/serialization/base_object.hpp>
#include <boost/shared_ptr.hpp>
#include "GeneralModelTransform.h"

namespace jiba
{
	/** \addtogroup inversion General routines for inversion */
	/* @{ */




	//! For the cross-gradient we sometimes have to extract to sections from the model vector that are not continuous in memory
	/*! For cases were we want to extract several sections from the generalized model parameters
	 * that are not continuous in memory, we cannot simply piece a transform
	 * from two different SectionTransform objects. For example for the cross-gradient between seismic and MT
	 * we currently store the model in the order slowness, density, conductivity. So we have to extract the
	 * first n model parameters and the last n model parameters and form a vector of length 2n that
	 * is fed to the cross-gradient function. This class generalizes this problem to an arbitrary number
	 * of sections. In principle each section can have a different length n1, n2, ... and the resulting
	 * transformed vector will have length n1+n2+...
	 */
	class MultiSectionTransform: public jiba::GeneralModelTransform
	{
	private:
		size_t length;
		std::vector<size_t> startindices;
		std::vector<size_t> endindices;
		std::vector<boost::shared_ptr<GeneralModelTransform> > Transforms;
		friend class boost::serialization::access;
		//! Provide serialization to be able to store objects and, more importantly for simpler MPI parallelization
		template<class Archive>
		void serialize(Archive & ar, const unsigned int version)
		{
			ar
					& boost::serialization::base_object<GeneralModelTransform>(
							*this);
			ar & length;
			ar & startindices;
			ar & endindices;
			ar & Transforms;
		}
		MultiSectionTransform() :
				length(0)
		{

		}
	public:

		//! We setup a clone function to have a virtual constructor and create polymorphic copies
		virtual MultiSectionTransform* clone() const
		{
			return new MultiSectionTransform(*this);
		}
		virtual jiba::rvec GeneralizedToPhysical(
				const jiba::rvec &FullModel) const
		{
			using boost::numeric::ublas::subrange;
			//we calculate the total length of the transformed vector from
			//the length of each section
			const size_t outlength = std::accumulate(endindices.begin(),
					endindices.end(), 0)
					- std::accumulate(startindices.begin(), startindices.end(),
							0);
			jiba::rvec Result(outlength);

			const size_t nsections = startindices.size();
			//the result will be continuous in memory and always
			//start at 0, each section will just follow after the other
			size_t resultstartindex = 0;
			for (size_t i = 0; i < nsections; ++i)
			{
				const size_t resultendindex = resultstartindex + endindices[i]
						- startindices[i];
				subrange(Result, resultstartindex, resultendindex) =
						Transforms[i]->GeneralizedToPhysical(
								subrange(FullModel, startindices[i],
										endindices[i]));
				resultstartindex = resultendindex;
			}
			return Result;
		}
		virtual jiba::rvec PhysicalToGeneralized(
				const jiba::rvec &FullModel) const
		{
			jiba::rvec Result(length);
			Result.clear();
			const size_t nsections = startindices.size();
			size_t currstart = 0;
			for (size_t i = 0; i < nsections; ++i)
			{
				subrange(Result, startindices[i], endindices[i]) =
						Transforms[i]->PhysicalToGeneralized(
								subrange(
										FullModel,
										currstart,
										currstart + endindices[i]
												- startindices[i]));
				currstart += endindices[i] - startindices[i];
			}
			return Result;
		}
		virtual jiba::rvec Derivative(const jiba::rvec &FullModel,
				const jiba::rvec &Derivative) const
		{
			using boost::numeric::ublas::subrange;
			jiba::rvec Result(length);
			Result.clear();
			const size_t nsections = startindices.size();
			size_t currstart = 0;
			for (size_t i = 0; i < nsections; ++i)
			{
				subrange(Result, startindices[i], endindices[i]) =
						Transforms[i]->Derivative(
								subrange(FullModel, startindices[i],
										endindices[i]),
								subrange(
										Derivative,
										currstart,
										currstart + endindices[i]
												- startindices[i]));
				currstart += endindices[i] - startindices[i];
			}
			return Result;
		}
		void AddSection(size_t startindex, size_t endindex,
				boost::shared_ptr<GeneralModelTransform> Trans)
		{
			startindices.push_back(startindex);
			endindices.push_back(endindex);
			Transforms.push_back(Trans);
		}
		MultiSectionTransform(size_t l) :
				length(l)
		{

		}
		MultiSectionTransform(size_t l, size_t startindex, size_t endindex,
				boost::shared_ptr<GeneralModelTransform> Trans) :
				length(l)
		{
			AddSection(startindex, endindex, Trans);
		}
		virtual ~MultiSectionTransform()
		{

		}
	};

/* @} */
}

#endif /* MULTISECTIONTRANSFORM_H_ */
