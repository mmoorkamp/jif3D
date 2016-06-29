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
    template<class ThreeDModelType>
    class J3DEXPORT MinMemGravMagCalculator: public jif3D::ThreeDGravMagCalculator<ThreeDModelType>
      {
    private:
      friend class access;
      //! Provide serialization to be able to store objects and, more importantly for simpler MPI parallelization
      template<class Archive>
      void serialize(Archive & ar, const unsigned int version)
        {
          ar
              & base_object<ThreeDGravMagCalculator<ThreeDModelType> >(
                  *this);
        }
    public:
      //! The implementation of the forward calculation
      virtual rvec Calculate(const ThreeDModelType &Model);
      //! We have to implement this function even though it does not do anything
      virtual void HandleSensitivities(const size_t measindex)
        {
        }
      //! The constructor takes a shared pointer to an implementation object
      MinMemGravMagCalculator(
          boost::shared_ptr<ThreeDGravMagImplementation<ThreeDModelType> > TheImp);
      virtual ~MinMemGravMagCalculator();
      };

    template<class ThreeDModelType>
    MinMemGravMagCalculator<ThreeDModelType>::MinMemGravMagCalculator(
        boost::shared_ptr<ThreeDGravMagImplementation<ThreeDModelType> > TheImp) :
        ThreeDGravMagCalculator<ThreeDModelType>(TheImp)
      {

      }

    template<class ThreeDModelType>
    MinMemGravMagCalculator<ThreeDModelType>::~MinMemGravMagCalculator()
      {

      }

    /*! For this simple case the calculator class only checks the model for consistency,
     * then calls the implementation object to do the calculation and returns the result
     * @param Model The Gravity model for the forward calculation
     * @return The resulting measurements (scalar or tensorial)
     */
    template<class ThreeDModelType>
    rvec MinMemGravMagCalculator<ThreeDModelType>::Calculate(const ThreeDModelType &Model)
      {
        ThreeDGravMagCalculator<ThreeDModelType>::SetCurrentSensitivities().resize(0, 0);
        return ThreeDGravMagCalculator<ThreeDModelType>::Calculate(Model);
      }
  /* @} */
  }

#endif /* MINMEMGRAVITYCALCULATOR_H_ */
