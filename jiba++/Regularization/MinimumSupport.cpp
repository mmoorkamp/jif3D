//============================================================================
// Name        : MatOpRegularization.cpp
// Author      : Jan 11, 2010
// Version     :
// Copyright   : 2010, mmoorkamp
//============================================================================

#include <boost/math/special_functions/pow.hpp>
#include "MinimumSupport.h"

namespace jiba
{
	namespace bm = boost::math;
	void MinimumSupport::ImplDataDifference(const jiba::rvec &Model,
			jiba::rvec &Diff)
	{
		RegFunc->CalcMisfit(Model);
		RegDiff = RegFunc->GetDataDifference();
		Diff.resize(Model.size(), false);
		const double b = beta * beta;
		const size_t nmod = Model.size();
		for (size_t i = 0; i < nmod; ++i)
		{
			double mag = bm::pow<2>(RegDiff(i)) + bm::pow<2>(RegDiff(nmod + i))
					+ bm::pow<2>(RegDiff(2 * nmod + i));
			//std::cout << RegDiff(i) << " " << RegDiff(2*i) << " " << RegDiff(3*i) << std::endl;
			Diff(i) = sqrt(mag / (b + mag));
		}
		//std::cout << Model.size() << " " << RegDiff.size() << std::endl;
	}

	jiba::rvec MinimumSupport::ImplGradient(const jiba::rvec &Model,
			const jiba::rvec &Diff)
	{
		jiba::rvec Grad = RegFunc->CalcGradient(Model);
		const double b = beta * beta;
		const size_t nmod = Model.size();
		for (size_t i = 0; i < Grad.size(); ++i)
		{
			std::cout << RegDiff(i) << " " << RegDiff(nmod + i) << " "
					<< RegDiff(2 * nmod + i) << std::endl;
			double mag = bm::pow<2>(RegDiff(i)) + bm::pow<2>(RegDiff(nmod + i))
					+ bm::pow<2>(RegDiff(2 * nmod + i));
			Grad(i) = 2 * b
					* (RegDiff(i) + RegDiff(nmod + i) + RegDiff(2 * nmod + i))
					/ bm::pow<2>(b + mag);
		}
		jiba::comp_mat Mat(
				RegFunc->GetXOperator() + RegFunc->GetYOperator()
						+ RegFunc->GetZOperator());
		return ublas::prod(trans(Mat), Grad);
	}

	MinimumSupport::MinimumSupport(
			boost::shared_ptr<jiba::MatOpRegularization> RF, double b) :
			beta(b), RegFunc(RF)
	{
	}

	MinimumSupport::~MinimumSupport()
	{
	}

} /* namespace jiba */
