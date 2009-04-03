//============================================================================
// Name        : FullSensitivityGravityCalculator.cpp
// Author      : Dec 12, 2008
// Version     :
// Copyright   : 2008, mmoorkamp
//============================================================================


#include "FullSensitivityGravityCalculator.h"

#ifdef HAVEATLAS
#include <boost/numeric/bindings/atlas/cblas2.hpp>
#include <boost/numeric/bindings/atlas/cblas3.hpp>
#endif

namespace jiba
  {

#ifdef HAVEATLAS
    namespace atlas = boost::numeric::bindings::atlas;
#endif

    FullSensitivityGravityCalculator::FullSensitivityGravityCalculator(
        boost::shared_ptr<ThreeDGravityImplementation> TheImp) :
      CachedGravityCalculator(TheImp)
      {

      }

    FullSensitivityGravityCalculator::~FullSensitivityGravityCalculator()
      {

      }

    void FullSensitivityGravityCalculator::HandleSensitivities(
        const size_t measindex)
      {
        const size_t startindex = measindex
            * Imp.get()->GetDataPerMeasurement();
        const size_t endindex = (measindex + 1)
            * Imp.get()->GetDataPerMeasurement();
        ublas::matrix_range<jiba::rmat> mr(Sensitivities, ublas::range(
            startindex, endindex), ublas::range(0, Sensitivities.size2()));
        mr = SetCurrentSensitivities();
      }

    rvec FullSensitivityGravityCalculator::CalculateNewModel(
        const ThreeDGravityModel &Model)
      {
        const size_t nmeas = Model.GetMeasPosX().size();
        const size_t ngrid = Model.GetDensities().num_elements();
        const size_t nmod = ngrid + Model.GetBackgroundThicknesses().size();
        Sensitivities.resize(nmeas * Imp.get()->GetDataPerMeasurement(), nmod,
            false);

        return Imp.get()->Calculate(Model, *this);
      }

    rvec FullSensitivityGravityCalculator::CalculateCachedResult(
        const ThreeDGravityModel &Model)
      {
        const size_t nmeas = Model.GetMeasPosX().size() * Imp.get()->GetDataPerMeasurement();
        const size_t ngrid = Model.GetDensities().num_elements();
        const size_t nmod = ngrid + Model.GetBackgroundThicknesses().size();

        assert(Sensitivities.size1() == nmeas);
        assert(Sensitivities.size2() == nmod);
        rvec result(nmeas );
        rvec DensVector(nmod);
        std::copy(Model.GetDensities().origin(), Model.GetDensities().origin()
            + ngrid, DensVector.begin());
        std::copy(Model.GetBackgroundDensities().begin(),
            Model.GetBackgroundDensities().end(), DensVector.begin() + ngrid);
#ifdef HAVEATLAS
        atlas::gemv(1.0, Sensitivities, DensVector, 0.0, result);
#else
        result = boost::numeric::ublas::prod(Sensitivities, DensVector);
#endif
        return result;
      }

    rvec FullSensitivityGravityCalculator::CachedLQDerivative(
        const ThreeDGravityModel &Model, const rvec &Misfit)
      {
        const size_t ngrid = Model.GetDensities().num_elements();
        const size_t nmod = ngrid + Model.GetBackgroundThicknesses().size();

        assert(Misfit.size() == Sensitivities.size1());
        assert(Sensitivities.size2() == nmod);
        rvec result(nmod);
#ifdef HAVEATLAS
        atlas::gemv(CblasTrans, 1.0, Sensitivities, Misfit, 0.0, result);
#else
        result = boost::numeric::ublas::prod(trans(Sensitivities),Misfit);
#endif
        return result;
      }
  }
