//============================================================================
// Name        : ThreeDGravityImplementation.h
// Author      : Dec 10, 2008
// Version     :
// Copyright   : 2008, mmoorkamp
//============================================================================


#ifndef THREEDGRAVITYIMPLEMENTATION_H_
#define THREEDGRAVITYIMPLEMENTATION_H_

#include "ThreeDGravityModel.h"
#include "../Global/VecMat.h"

namespace jiba
  {
    /** \addtogroup gravity Gravity forward modeling, display and inversion */
    /* @{ */
    //this is just a forward declaration to avoid circular inclusions
    class ThreeDGravityCalculator;
    //! The base class that provides the interface for the numerical implementation of the gravity forward calculations.
    /*! The calculation of the forward response is split into two class hierarchies that
     * have to be used in conjunction. The classes derived from ThreeDGravityImplementation
     * are responsible for the numerical implementation details of the calculation, i.e. using
     * parallelization with openmp or performing calculations on a graphics card with cuda, and
     * whether FTG or scalar data is calculated. These classes are not directly visible to the
     * user, but only through the calculator classes that determine whether sensitivity information
     * is stored, whether calculations are cached with stored sensitivity information and how this
     * information is processed and stored.
     *
     * This type of design resembles the bridge pattern and allows to freely combine optimized implementations
     * for different platforms with different sensitivity handlings.
     */
    class ThreeDGravityImplementation
      {
    private:
      void CacheGeometry(const ThreeDGravityModel &Model);
      //! Calculate the response of the 1D background for a single measurement, this function has to be implemented in the derived class.
      virtual rvec CalcBackground(const size_t measindex, const double xwidth, const double ywidth,
          const double zwidth, const ThreeDGravityModel &Model,
          rmat &Sensitivities) = 0;
      //! Calculate the response of the gridded domain for a single measurement, this function has to be implemented in the derived class.
      virtual rvec CalcGridded(const size_t measindex, const ThreeDGravityModel &Model,
          rmat &Sensitivities) = 0;
    protected:
      // The access functions for the coordinates are not thread safe
      // so we cache the values once in the default implementation of calculate.
      // Because CalcBackground and CalcGridded can only be called from within
      // there we are guaranteed correct values in the implementation of those functions
      ThreeDGravityModel::t3DModelDim XCoord;
      ThreeDGravityModel::t3DModelDim YCoord;
      ThreeDGravityModel::t3DModelDim ZCoord;
      ThreeDGravityModel::t3DModelDim XSizes;
      ThreeDGravityModel::t3DModelDim YSizes;
      ThreeDGravityModel::t3DModelDim ZSizes;
    public:
      //! We can implement tensor and scalar calculations in the derived classes, this function returns how many data values a single measurement yields
      virtual size_t GetDataPerMeasurement() = 0;
      //! For a given Model calculate the forward response for all measurements and return it as a real vector, the calculator object is passed to process the sensitivity information
      virtual rvec Calculate(const ThreeDGravityModel &Model,
          ThreeDGravityCalculator &Calculator);
      virtual rvec LQDerivative(const ThreeDGravityModel &Model, const rvec &Misfit,ThreeDGravityCalculator &Calculator);
      ThreeDGravityImplementation();
      virtual ~ThreeDGravityImplementation();
      };
  /* @} */
  }
#endif /* THREEDGRAVITYIMPLEMENTATION_H_ */
