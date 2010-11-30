//============================================================================
// Name        : ThreeDGravityCalculator.h
// Author      : Dec 10, 2008
// Version     :
// Copyright   : 2008, mmoorkamp
//============================================================================


#ifndef THREEDGRAVITYCALCULATOR_H_
#define THREEDGRAVITYCALCULATOR_H_

#include "ThreeDGravityModel.h"
#include "ThreeDGravityImplementation.h"
#include <boost/shared_ptr.hpp>

namespace jiba
  {
    /** \addtogroup gravity Gravity forward modeling, display and inversion */
    /* @{ */
    //! The base class for all calculator objects, these handle how the sensitivities are stored and processed
    /*! The ThreeDGravityCalculator class is the base class that provides
     * the user interface to the gravity forward calculation. It uses a
     * ThreeDGravityImplementation object to do the actual forward calculation.
     * The derived classes determine how the sensitivity information
     * is handled, whether it is discarded to save memory, fully stored, or
     * compressed in some way. Also they decide in how far calculations are
     * cached using sensitivity information.
     *
     * Within the forward calculation for each measurement the implementation class
     * calls the member function HandleSensitivities() with the index for the measurement
     * to give the calculator class the possibility to store, process or discard the
     * current row of the sensitivity matrix. This avoids having to store the complete matrix
     * in cases where some compression is applied.
     *
     * Note also that the calculator class hierarchy is independent of whether
     * scalar or FTG data or a derived quantity is calculated. This is purely determined by the Implementation object.
     * Therefore the data is returned as a single vector. If we calculate scalar data it has one entry
     * per measurement point in the model. For FTG data it contains 9 consecutive entries with the matrix
     * elements per measurement.
     *
     * We can set a transformation to directly calculate derived quantities and the associated gradient. This transform
     * is forwarded to the implementation object and all returned data and gradients will be with respect to this transformation.
     */
    class ThreeDGravityCalculator
      {
    public:
      typedef ThreeDGravityModel ModelType;
    private:
      /*! We need a structure to hold the sensitivities for
       * the current measurement that can be passed to
       * the implementation object. The derived classes
       * decide how to handle this information
       */
      rmat CurrentSensitivities;
    protected:
      //! In some cases we might want to apply a transformation to the data, e.g. FTG to an invariant
      boost::shared_ptr<VectorTransform> Transform;
      //! The shared pointer to the implementation object that does the actual calculation
      boost::shared_ptr<ThreeDGravityImplementation> Imp;
      //! Check the the information in the model is consistent, i.e. corresponding vectors have the same size
      void CheckModelConsistency(const ModelType &Model);
    public:
      //! Assign an object that performs a transformation on the data, e.g. FTG tensor to an invariant
      void SetDataTransform(boost::shared_ptr<VectorTransform> DataTransform)
        {
          Transform = DataTransform;
          Imp->SetDataTransform(DataTransform);
        }
      //! Calculate the forward response of the given model, this simple implementation just forwards the call to the implementation class
      virtual rvec Calculate(const ModelType &Model);
      //! Get the least squares derivative \f$ \partial O/ \partial \f$ of a least squares objective function \f$ O = \sum (d^{obs} - d^{pred})^2 \f$
      virtual rvec LQDerivative(const ModelType &Model, const rvec &Misfit);
      //! Read and write access to the sensitivity information for the current measurement, only intended for implementation classes
      rmat &SetCurrentSensitivities()
        {
          return CurrentSensitivities;
        }
      //! In some cases we need to know the amount of data we get per measurement, this information is stored in the implementation class. It is 1 for scalar and at the moment 9 for FTG data
      size_t GetDataPerMeasurement()
        {
          return Imp->GetDataPerMeasurement();
        }
      //! Process the sensitivity information for the current measurement, called from the implementation class
      virtual void HandleSensitivities(const size_t measindex) = 0;
      //! This class is useless without an implementation object so we have to pass one to the constructor
      ThreeDGravityCalculator(
          boost::shared_ptr<ThreeDGravityImplementation> TheImp);
      virtual ~ThreeDGravityCalculator();
      };

  /* @} */
  }

#endif /* THREEDGRAVITYCALCULATOR_H_ */
