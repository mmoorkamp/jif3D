//============================================================================
// Name        : ThreeDGravityImplementation.h
// Author      : Dec 10, 2008
// Version     :
// Copyright   : 2008, mmoorkamp
//============================================================================


#ifndef THREEDGRAVITYIMPLEMENTATION_H_
#define THREEDGRAVITYIMPLEMENTATION_H_

#include <boost/shared_ptr.hpp>
#include "ThreeDGravityModel.h"
#include "../Global/VecMat.h"
#include "../Global/VectorTransform.h"

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
      boost::shared_ptr<VectorTransform> Transform;
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
      //! The cached values of the x-coordinates in m for each cell
      ThreeDGravityModel::t3DModelDim XCoord;
      //! The cached values of the y-coordinates in m for each cell
      ThreeDGravityModel::t3DModelDim YCoord;
      //! The cached values of the z-coordinates in m for each cell
      ThreeDGravityModel::t3DModelDim ZCoord;
      //! The cached values for the size of each cell in x-direction
      ThreeDGravityModel::t3DModelDim XSizes;
      //! The cached values for the size of each cell in y-direction
      ThreeDGravityModel::t3DModelDim YSizes;
      //! The cached values for the size of each cell in z-direction
      ThreeDGravityModel::t3DModelDim ZSizes;
    public:
      //! Returns the number of data before any transformation is applied
      virtual size_t RawDataPerMeasurement() = 0;
      //! Set a transformation class that should be applied to any data and gradient in the calculation
      void SetDataTransform(boost::shared_ptr<VectorTransform> DataTransform)
        {
          Transform = DataTransform;
        }
      //! We can implement tensor and scalar calculations in the derived classes, this function returns how many data values a single measurement yields and considers any transformation
      size_t GetDataPerMeasurement()
        {
          return Transform ? Transform->GetOutputSize() : RawDataPerMeasurement();
        }
      //! For a given Model calculate the forward response for all measurements and return it as a real vector, the calculator object is passed to process the sensitivity information
      virtual rvec Calculate(const ThreeDGravityModel &Model,
          ThreeDGravityCalculator &Calculator);
      //! Calculate the least-squres derivative vector for the given model and vector
      virtual rvec LQDerivative(const ThreeDGravityModel &Model, const rvec &Misfit);
      ThreeDGravityImplementation();
      virtual ~ThreeDGravityImplementation();
      };
  /* @} */
  }
#endif /* THREEDGRAVITYIMPLEMENTATION_H_ */
