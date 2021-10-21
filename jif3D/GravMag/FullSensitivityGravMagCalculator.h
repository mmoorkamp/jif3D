//============================================================================
// Name        : FullSensitivityGravityCalculator.h
// Author      : Dec 12, 2008
// Version     :
// Copyright   : 2008, mmoorkamp
//============================================================================

#ifndef FULLSENSITIVITYGRAVITYCALCULATOR_H_
#define FULLSENSITIVITYGRAVITYCALCULATOR_H_

#include "../Global/Serialization.h"
#include "../Global/Jif3DGlobal.h"
#include "../Global/FatalException.h"
#include "CachedGravMagCalculator.h"

#include <iostream>
namespace jif3D
  {
    /** \addtogroup gravity Gravity forward modeling, display and inversion */
    /* @{ */
    //! A calculator class that stores the full sensitivity matrix for inversion or fast consecutive forward calculations
    /*! This class stores the complete sensitivity matrix with all entries.
     * This allows to perform direct linear inversion and makes series of forward calculations
     * with identical geometries but changing density values extremely fast.
     * The disadvantage of this approach is that large models take up huge amounts of memory.
     * The actual management of the caching of calculations is performed by the base class
     * CachedGravityCalculator, this class only manages the storage and how cached results
     * are calculated.
     *
     * Important note: The stored sensitivities are the raw sensitivities irrespective of any
     * transformation. This is necessary to calculate correct gradients for transformed data.
     */
    template<class PotentialDataType>
    class J3DEXPORT FullSensitivityGravMagCalculator: public jif3D::CachedGravMagCalculator<
        PotentialDataType>
      {
    private:
      typedef typename PotentialDataType::ModelType ThreeDModelType;

      //! The full sensitivity matrix G for the gravity measurements, we use it to calculate the forward response through d = G*m
      rmat Sensitivities;
      //! Calculates the raw gravity/FTG data from cached sensitivities without applying any transformation
      virtual rvec CalculateRawData(const ThreeDModelType &Model,
          const PotentialDataType &Data);
      //! Calculates the raw gravity/FTG derivative from cached sensitivities without applying any transformation
      virtual rvec CalculateRawLQDerivative(const ThreeDModelType &Model,
          const PotentialDataType &Data, const rvec &Misfit);
      //! Calculate a new model when no cached information is available or valid
      virtual rvec CalculateNewModel(const ThreeDModelType &Model,
          const PotentialDataType &Data) override;
      //! calculate Data with applied transformation when cached sensitivities are still valid
      virtual rvec CalculateCachedResult(const ThreeDModelType &Model,
          const PotentialDataType &Data) override;
      //! calculate derivative of least squares objective function with applied transformation when cached sensitivities are still valid
      virtual rvec CachedLQDerivative(const ThreeDModelType &Model,
          const PotentialDataType &Data, const rvec &Misfit) override;
      friend class access;
      //! Provide serialization to be able to store objects and, more importantly for simpler MPI parallelization
      template<class Archive>
      void serialize(Archive & ar, const unsigned int version)
        {
          ar & base_object<ThreeDGravMagCalculator>(*this);
          ar & Sensitivities;
        }
    public:
      //! return a read only copy of the sensitivity matrix, this guarantees that cache information is preserved
      const rmat &GetSensitivities() const
        {
          return Sensitivities;
        }
      //! For efficiency we sometimes operate directly on the sensitivities, as we don't have guaranteed cache information, we enforce recalculation
      rmat &SetSensitivities()
        {
          CachedGravMagCalculator<PotentialDataType>::InvalidateCache();
          return Sensitivities;
        }
      //! This function is called by the implementation classes and allows to integrate the results from a single measurement
      virtual void HandleSensitivities(const size_t measindex) override;
      //! The constructor takes a shared pointer to an implementation object
      FullSensitivityGravMagCalculator(
          boost::shared_ptr<ThreeDGravMagImplementation<PotentialDataType> > TheImp);
      virtual ~FullSensitivityGravMagCalculator();
      };
    template<class PotentialDataType>
    FullSensitivityGravMagCalculator<PotentialDataType>::FullSensitivityGravMagCalculator(
        boost::shared_ptr<ThreeDGravMagImplementation<PotentialDataType> > TheImp) :
        CachedGravMagCalculator<PotentialDataType>(TheImp), Sensitivities()
      {
std::cout << "Storing sensitivities in memory " << std::endl;
      }

    template<class PotentialDataType>
    FullSensitivityGravMagCalculator<PotentialDataType>::~FullSensitivityGravMagCalculator()
      {

      }

    template<class PotentialDataType>
    void FullSensitivityGravMagCalculator<PotentialDataType>::HandleSensitivities(
        const size_t measindex)
      {
        //we have to identify the correct rows in the full sensitivity
        //matrix where we want to store the current sensitivity information
        const size_t startindex =
            measindex
                * ThreeDGravMagCalculator<PotentialDataType>::Imp.get()->RawDataPerMeasurement();
        const size_t endindex =
            (measindex + 1)
                * ThreeDGravMagCalculator<PotentialDataType>::Imp.get()->RawDataPerMeasurement();
        //construct a range that we can assign the current sensitivities to
        ublas::matrix_range<jif3D::rmat> mr(Sensitivities,
            ublas::range(startindex, endindex), ublas::range(0, Sensitivities.size2()));
        mr = ThreeDGravMagCalculator<PotentialDataType>::SetCurrentSensitivities();
      }

    template<class PotentialDataType>
    rvec FullSensitivityGravMagCalculator<PotentialDataType>::CalculateNewModel(
        const ThreeDModelType &Model, const PotentialDataType &Data)
      {
        //allocate enough memory for the sensitivities
        const size_t nmeas = Data.GetMeasPosX().size();
        const size_t nmod = Model.GetNModelParm();
        Sensitivities.resize(
            nmeas
                * ThreeDGravMagCalculator<PotentialDataType>::Imp.get()->RawDataPerMeasurement(),
            nmod, false);
        //then forward the call to the implementation object
        return ThreeDGravMagCalculator<PotentialDataType>::Imp.get()->Calculate(Model,
            Data, *this);
      }

    //When we have sensitivity information we first
    //calculate the untransformed data and then
    //apply the appropriate transformations for
    //the gradient or the transformed data
    template<class PotentialDataType>
    rvec FullSensitivityGravMagCalculator<PotentialDataType>::CalculateRawData(
        const ThreeDModelType &Model, const PotentialDataType &Data)
      {
        const size_t nmeas =
            Data.GetMeasPosX().size()
                * ThreeDGravMagCalculator<PotentialDataType>::Imp.get()->RawDataPerMeasurement();
        const size_t nmod = Model.GetNModelParm();
        if (Sensitivities.size1() != nmeas || Sensitivities.size2() != nmod)
          {
            throw jif3D::FatalException(
                "Size of sensitivity matrix does not match model configuration ! ",
                __FILE__, __LINE__);
          }
        //perform the vector matrix product to get the raw data
        return boost::numeric::ublas::prec_prod(Sensitivities, Model.GetModelParameters());
      }

    template<class PotentialDataType>
    rvec FullSensitivityGravMagCalculator<PotentialDataType>::CalculateCachedResult(
        const ThreeDModelType &Model, const PotentialDataType &Data)
      {
        //if we have to transform the data
        if (ThreeDGravMagCalculator<PotentialDataType>::Transform)
          {
            //we apply the transform to the calculated raw data
            return ApplyTransform(CalculateRawData(Model, Data),
                *ThreeDGravMagCalculator<PotentialDataType>::Transform);
          }
        //otherwise we just return the raw data
        return CalculateRawData(Model, Data);
      }

    template<class PotentialDataType>
    rvec FullSensitivityGravMagCalculator<PotentialDataType>::CachedLQDerivative(
        const ThreeDModelType &Model, const PotentialDataType &Data, const rvec &Misfit)
      {
        //first we calculate the raw data, the transformation might depend on this data
        rvec Result(CalculateRawData(Model, Data));
        rvec ProcessedMisfit;
        //if we don't have to apply a transform
        if (!ThreeDGravMagCalculator<PotentialDataType>::Transform)
          {
            // we just copy the misfit for further processing
            ProcessedMisfit = Misfit;
          }
        else
          {
            //otherwise we have to process each segment
            //of the misift of the transformed data
            //to calculate the correct gradient
            const size_t nin =
                ThreeDGravMagCalculator<PotentialDataType>::Transform->GetInputSize();
            const size_t nout =
                ThreeDGravMagCalculator<PotentialDataType>::Transform->GetOutputSize();
            //the sensitivities are for the raw data
            //so our new processed misfit has to have the same size as the raw data
            ProcessedMisfit.resize(Misfit.size() / nout * nin);
            for (size_t i = 0; i < Misfit.size(); i += nout)
              {
                size_t outindex = i / nout * nin;
                ublas::vector_range<jif3D::rvec> OutRange(ProcessedMisfit,
                    ublas::range(outindex, outindex + nin));
                ublas::vector_range<const jif3D::rvec> InRange(Misfit,
                    ublas::range(i, i + nout));
                ublas::vector_range<const jif3D::rvec> DataRange(Result,
                    ublas::range(outindex, outindex + nin));
                OutRange = ublas::prod(
                    trans(
                        ThreeDGravMagCalculator<PotentialDataType>::Transform->Derivative(
                            DataRange)), InRange);
              }

          }

        return CalculateRawLQDerivative(Model, Data, ProcessedMisfit);
      }

    template<class PotentialDataType>
    rvec FullSensitivityGravMagCalculator<PotentialDataType>::CalculateRawLQDerivative(
        const ThreeDModelType &Model, const PotentialDataType &Data, const rvec &Misfit)
      {
        //when we have the sensitivities, the derivative of the objective function
        //is simply J^T * delta d
        return 2.0 * boost::numeric::ublas::prod(trans(Sensitivities), Misfit);
      }
  /* @} */
  }

#endif /* FULLSENSITIVITYGRAVITYCALCULATOR_H_ */
