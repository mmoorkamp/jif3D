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
      CachedGravityCalculator(TheImp), Sensitivities()
      {

      }

    FullSensitivityGravityCalculator::~FullSensitivityGravityCalculator()
      {

      }

    void FullSensitivityGravityCalculator::HandleSensitivities(
        const size_t measindex)
      {
        //we have to identify the correct rows in the full sensitivity
        //matrix where we want to store the current sensitivity information
        const size_t startindex = measindex
            * Imp.get()->RawDataPerMeasurement();
        const size_t endindex = (measindex + 1)
            * Imp.get()->RawDataPerMeasurement();
        //construct a range that we can assign the current sensitivities to
        ublas::matrix_range<jiba::rmat> mr(Sensitivities, ublas::range(
            startindex, endindex), ublas::range(0, Sensitivities.size2()));
        mr = SetCurrentSensitivities();
      }

    rvec FullSensitivityGravityCalculator::CalculateNewModel(
        const ThreeDGravityModel &Model)
      {
        //allocate enough memory for the sensitivities
        const size_t nmeas = Model.GetMeasPosX().size();
        const size_t ngrid = Model.GetDensities().num_elements();
        const size_t nmod = ngrid + Model.GetBackgroundThicknesses().size();
        Sensitivities.resize(nmeas * Imp.get()->RawDataPerMeasurement(), nmod,
            false);
        //then forward the call to the implementation object
        return Imp.get()->Calculate(Model, *this);
      }

    //When we have sensitivity information we first
    //calculate the untransformed data and then
    //apply the appropriate transformations for
    //the gradient or the transformed data
    rvec FullSensitivityGravityCalculator::CalculateRawData(
        const ThreeDGravityModel &Model)
      {
        const size_t nmeas = Model.GetMeasPosX().size()
            * Imp.get()->RawDataPerMeasurement();
        const size_t ngrid = Model.GetDensities().num_elements();
        const size_t nmod = ngrid + Model.GetBackgroundThicknesses().size();

        assert(Sensitivities.size1() == nmeas);
        assert(Sensitivities.size2() == nmod);
        rvec DensVector(nmod);
        //copy the 3D model structure and the background into a vector of densities
        std::copy(Model.GetDensities().origin(), Model.GetDensities().origin()
            + ngrid, DensVector.begin());
        std::copy(Model.GetBackgroundDensities().begin(),
            Model.GetBackgroundDensities().end(), DensVector.begin() + ngrid);
        //perform the vector matrix product to get the raw data
        //depending on whether atlas is available on the system
        //we use the fast atlas version or the slower native ublas code
#ifdef HAVEATLAS
        rvec result(nmeas);
        atlas::gemv(1.0, Sensitivities, DensVector, 0.0, result);
        return result;
#else
        return boost::numeric::ublas::prod(Sensitivities, DensVector);
#endif

      }

    rvec FullSensitivityGravityCalculator::CalculateCachedResult(
        const ThreeDGravityModel &Model)
      {
        //if we have to transform the data
        if (Transform)
          {
            //we apply the transform to the calculated raw data
            return ApplyTransform(CalculateRawData(Model), *Transform);
          }
        //otherwise we just return the raw data
        return CalculateRawData(Model);
      }

    rvec FullSensitivityGravityCalculator::CachedLQDerivative(
        const ThreeDGravityModel &Model, const rvec &Misfit)
      {
        //first we calculate the raw data, the transformation might depend on this data
        rvec Data(CalculateRawData(Model));
        rvec ProcessedMisfit;
        //if we don't have to apply a transform
        if (!Transform)
          {
            // we just copy the misfit for further processing
            ProcessedMisfit = Misfit;
          }
        else
          {
            //otherwise we have to process each segment
            //of the misift of the transformed data
            //to calculate the correct gradient
            const size_t nin = Transform->GetInputSize();
            const size_t nout = Transform->GetOutputSize();
            //the sensitivities are for the raw data
            //so our new processed misfit has to have the same size as the raw data
            ProcessedMisfit.resize(Misfit.size() / nout * nin);
            for (size_t i = 0; i < Misfit.size(); i += nout)
              {
                size_t outindex = i / nout * nin;
                ublas::vector_range<jiba::rvec> OutRange(ProcessedMisfit,
                    ublas::range(outindex, outindex + nin));
                ublas::vector_range<const jiba::rvec> InRange(Misfit,
                    ublas::range(i, i + nout));
                ublas::vector_range<const jiba::rvec> DataRange(Data,
                    ublas::range(outindex, outindex + nin));
                OutRange = ublas::prod(trans(Transform->Derivative(DataRange)),
                    InRange);
              }

          }

        return CalculateRawLQDerivative(Model, ProcessedMisfit);
      }

    rvec FullSensitivityGravityCalculator::CalculateRawLQDerivative(
        const ThreeDGravityModel &Model, const rvec &Misfit)
      {
        //when we have the sensitivities, the derivative of the objective function
        //is simply J^T * delta d, as above we use the fast atlas matrix
        //vector operation when available and the slower ublas version otherwise
        assert(Misfit.size() == Sensitivities.size1());
#ifdef HAVEATLAS
        rvec result(Misfit.size());
        atlas::gemv(CblasTrans, 2.0, Sensitivities, Misfit, 0.0, result);
        return result;
#else
        return 2.0 * boost::numeric::ublas::prod(trans(Sensitivities), Misfit);
#endif
      }
  }
