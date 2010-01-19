//============================================================================
// Name        : TomographyObjective.h
// Author      : May 5, 2009
// Version     :
// Copyright   : 2009, mmoorkamp
//============================================================================


#ifndef TOMOGRAPHYOBJECTIVE_H_
#define TOMOGRAPHYOBJECTIVE_H_

#include "../ModelBase/ModelRefiner.h"
#include "../Inversion/ObjectiveFunction.h"
#include "TomographyCalculator.h"
#include "ThreeDSeismicModel.h"

namespace jiba
  {
    /** \addtogroup tomo Seismic tomography classes and functions */
    /* @{ */

    //! The objective function for seismic tomography
    /*! This class implements the objective function for 3D seismic refraction calculations.
     * We use a TomographyCalculator object to calculate the travel times for each measurement
     * and the the associated rays. From the rays we construct the sensitivities and the
     * gradient of the objective function with respect to the model parameters.
     *
     */
    class TomographyObjective: public ObjectiveFunction
      {
    private:
      // The slowness model we use for the forward calculations, keeps track of model and measurement geometry
      jiba::ThreeDSeismicModel FineSlownessModel;
      jiba::ThreeDSeismicModel CoarseSlownessModel;
      //We might use a different discretization for forward modeling and inversion
      jiba::ModelRefiner Refiner;
      // The observed travel times
      jiba::rvec ObservedData;
      // The calculator object
      TomographyCalculator Calculator;
      //The implementation of the function to calculate the difference between observed and synthetic data
      virtual void
      ImplDataDifference(const jiba::rvec &Model, jiba::rvec &Diff);
      //The implementation of the gradient calculation
      virtual jiba::rvec ImplGradient(const jiba::rvec &Model,
          const jiba::rvec &Diff);
      //! So far transformations have no effect for travel time data !
      virtual void SetDataTransformAction()
        {
          //TODO Implement transformation if necessary
        }
    public:
      //! Set the observed travel times
      void SetObservedData(const jiba::rvec &Data)
        {
          ObservedData = Data;
        }
      //! Set the model that contains the grid geometry and measurement and shot positions
      void SetFineModelGeometry(const jiba::ThreeDSeismicModel &Model)
        {
          FineSlownessModel = Model;
        }
      void SetCoarseModelGeometry(const jiba::ThreeDSeismicModel &Model)
        {
          CoarseSlownessModel = Model;
        }
      TomographyObjective();
      virtual ~TomographyObjective();
      };
  /* @} */
  }

#endif /* TOMOGRAPHYOBJECTIVE_H_ */
