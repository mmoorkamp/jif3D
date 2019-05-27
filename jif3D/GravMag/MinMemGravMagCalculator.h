//============================================================================
// Name        : MinMemGravityCalculator.h
// Author      : Dec 10, 2008
// Version     :
// Copyright   : 2008, mmoorkamp
//============================================================================

#ifndef MINMEMGRAVITYCALCULATOR_H_
#define MINMEMGRAVITYCALCULATOR_H_

#include "../Global/Serialization.h"
#include "../Global/Jif3DGlobal.h"
#include "ThreeDGravMagCalculator.h"
#include <boost/shared_ptr.hpp>

namespace jif3D
  {
    /** \addtogroup gravity Gravity forward modeling, display and inversion */
    /* @{ */
    //! A calculator class that uses the minimum amount of memory, no caching, no sensitivity information
    /*! This is a calculator class that uses only the minimum amount of memory. It is suitable
     * when the forward response for a certain model geometry only has to be calculated once
     * or if the model is so big that the sensitivity matrix cannot be stored in memory or on disk any more.
     */
    template<class PotentialDataType>
    class J3DEXPORT MinMemGravMagCalculator: public jif3D::ThreeDGravMagCalculator<
        PotentialDataType>
      {
    private:
      typedef typename PotentialDataType::ModelType ThreeDModelType;

      friend class access;
      //! Provide serialization to be able to store objects and, more importantly for simpler MPI parallelization
      template<class Archive>
      void serialize(Archive & ar, const unsigned int version)
        {
          ar & base_object<ThreeDGravMagCalculator<PotentialDataType> >(*this);
        }
    public:
      //! The implementation of the forward calculation
      virtual rvec Calculate(const ThreeDModelType &Model, const PotentialDataType &Data)
          override;
      //! We have to implement this function even though it does not do anything
      virtual void HandleSensitivities(const size_t measindex) override
        {
        }
      //! The constructor takes a shared pointer to an implementation object
      MinMemGravMagCalculator(
          boost::shared_ptr<ThreeDGravMagImplementation<PotentialDataType> > TheImp);
      virtual ~MinMemGravMagCalculator();
      };

    template<class PotentialDataType>
    MinMemGravMagCalculator<PotentialDataType>::MinMemGravMagCalculator(
        boost::shared_ptr<ThreeDGravMagImplementation<PotentialDataType> > TheImp) :
        ThreeDGravMagCalculator<PotentialDataType>(TheImp)
      {

      }

    template<class PotentialDataType>
    MinMemGravMagCalculator<PotentialDataType>::~MinMemGravMagCalculator()
      {

      }

    /*! For this simple case the calculator class only checks the model for consistency,
     * then calls the implementation object to do the calculation and returns the result
     * @param Model The Gravity model for the forward calculation
     * @return The resulting measurements (scalar or tensorial)
     */
    template<class PotentialDataType>
    rvec MinMemGravMagCalculator<PotentialDataType>::Calculate(
        const ThreeDModelType &Model, const PotentialDataType &Data)
      {
        ThreeDGravMagCalculator<PotentialDataType>::SetCurrentSensitivities().resize(0,
            0);
        return ThreeDGravMagCalculator<PotentialDataType>::Calculate(Model, Data);
      }
  /* @} */
  }

#endif /* MINMEMGRAVITYCALCULATOR_H_ */
