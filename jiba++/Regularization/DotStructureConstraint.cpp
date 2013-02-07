/*
 * DotStructureConstraint.cpp
 *
 *  Created on: Feb 5, 2013
 *      Author: mmoorkamp
 */
#include <boost/math/special_functions/pow.hpp>
#include "DotStructureConstraint.h"
#include "../Global/NumUtil.h"

namespace jiba
{
	namespace bm = boost::math;
	void DotStructureConstraint::ImplDataDifference(const jiba::rvec &Model,
			jiba::rvec &Diff)
	{
		const size_t nmod = Model.size();
		//We have two complete models and three components of the cross gradient
		//we calculate the cross-gradient for each cell of the grid
		const size_t ndiff = nmod / 2;
		Diff.resize(ndiff);
		//we need the spatial gradients of both models, so we use the
		//appropriate objective function objects for this
		FirstGradient.CalcMisfit(
				ublas::vector_range<const jiba::rvec>(Model,
						ublas::range(0, ndiff)));
		SecondGradient.CalcMisfit(
				ublas::vector_range<const jiba::rvec>(Model,
						ublas::range(ndiff, nmod)));
		for (size_t i = 0; i < ndiff; ++i)
		{
			Diff(i) =
					(bm::pow<2>(FirstGradient.GetDataDifference()(i))
							+ bm::pow<2>(
									FirstGradient.GetDataDifference()(
											i + ndiff))
							+ bm::pow<2>(
									FirstGradient.GetDataDifference()(i + nmod)))
							* (bm::pow<2>(SecondGradient.GetDataDifference()(i))
									+ bm::pow<2>(
											SecondGradient.GetDataDifference()(
													i + ndiff))
									+ bm::pow<2>(
											SecondGradient.GetDataDifference()(
													i + nmod)))
							- bm::pow<2>(
									FirstGradient.GetDataDifference()(i)
											* SecondGradient.GetDataDifference()(
													i)
											+ FirstGradient.GetDataDifference()(
													i + ndiff)
													* SecondGradient.GetDataDifference()(
															i + ndiff)
											+ FirstGradient.GetDataDifference()(
													i + nmod)
													* SecondGradient.GetDataDifference()(
															i + nmod));
			Diff(i) = sign(Diff(i)) * sqrt(std::abs(Diff(i)));

		}
	}

	jiba::rvec DotStructureConstraint::ImplGradient(const jiba::rvec &Model,
			const jiba::rvec &Diff)
	{
		const size_t nmod = Model.size();
		const size_t halfmod = nmod / 2;
		FirstGradient.CalcGradient(
				ublas::vector_range<const jiba::rvec>(Model,
						ublas::range(0, halfmod)));
		SecondGradient.CalcGradient(
				ublas::vector_range<const jiba::rvec>(Model,
						ublas::range(halfmod, nmod)));
		jiba::rvec Gradient(nmod);
		Gradient.clear();
		for (size_t i = 0; i < halfmod; ++i)
		{
			Gradient(i) =
					(FirstGradient.GetXGrad()(i)
							* FirstGradient.GetDataDifference()(i)
							+ FirstGradient.GetYGrad()(i)
									* FirstGradient.GetDataDifference()(
											i + halfmod)
							+ FirstGradient.GetZGrad()(i)
							* FirstGradient.GetDataDifference()(i + nmod))
							* (bm::pow<2>(SecondGradient.GetDataDifference()(i))
									+ bm::pow<2>(
											SecondGradient.GetDataDifference()(
													i + halfmod))
									+ bm::pow<2>(
											SecondGradient.GetDataDifference()(
													i + nmod)))
							- (FirstGradient.GetDataDifference()(i)
									* SecondGradient.GetDataDifference()(i)
									+ FirstGradient.GetDataDifference()(
											i + halfmod)
											* SecondGradient.GetDataDifference()(
													i + halfmod)
									+ FirstGradient.GetDataDifference()(
											i + nmod)
											* SecondGradient.GetDataDifference()(
													i + nmod))

									* (FirstGradient.GetXGrad()(i)
											* SecondGradient.GetDataDifference()(
													i)
											+ FirstGradient.GetYGrad()(i)
													* SecondGradient.GetDataDifference()(
															i + halfmod)
											+ FirstGradient.GetZGrad()(i)
													* SecondGradient.GetDataDifference()(
															i + nmod));

			Gradient(i + halfmod) =
					(SecondGradient.GetXGrad()(i)
							* SecondGradient.GetDataDifference()(i)
							+ SecondGradient.GetYGrad()(i)
									* SecondGradient.GetDataDifference()(
											i + halfmod)
							+ SecondGradient.GetZGrad()(i)
									* SecondGradient.GetDataDifference()(
											i + nmod))
							* (bm::pow<2>(FirstGradient.GetDataDifference()(i))
									+ bm::pow<2>(
											FirstGradient.GetDataDifference()(
													i + halfmod))
									+ bm::pow<2>(
											FirstGradient.GetDataDifference()(
													i + nmod)))
							- (SecondGradient.GetDataDifference()(i)
									* FirstGradient.GetDataDifference()(i)
									+ SecondGradient.GetDataDifference()(
											i + halfmod)
											* FirstGradient.GetDataDifference()(
													i + halfmod)
									+ SecondGradient.GetDataDifference()(
											i + nmod)
											* FirstGradient.GetDataDifference()(
													i + nmod))

									* (SecondGradient.GetXGrad()(i)
											* FirstGradient.GetDataDifference()(
													i)
											+ SecondGradient.GetYGrad()(i)
													* FirstGradient.GetDataDifference()(
															i + halfmod)
											+ SecondGradient.GetZGrad()(i)
													* FirstGradient.GetDataDifference()(
															i + nmod));
		}
		return 2.0 * Gradient;
	}

	DotStructureConstraint::~DotStructureConstraint()
	{
		// TODO Auto-generated destructor stub
	}

} /* namespace jiba */
