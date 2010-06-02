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
      //! The slowness model we use for the forward calculations, keeps track of model and measurement geometry
      jiba::ThreeDSeismicModel FineSlownessModel;
      //! The goemetry of the slowness model used in the inversion
      jiba::ThreeDSeismicModel CoarseSlownessModel;
      //! We might use a different discretization for forward modeling and inversion
      /*! To ensure numerically correct travel times we might need a model discretization
       * that is much finer than the inversion grid we want to use. The Refiner object
       * takes CoarseSlownessModel which is the model with a geometry as intended in the
       * inversion and refines it using the geometry of FineSlownessModel. This model
       * is then used for the forward calculation.
       */
      jiba::ModelRefiner Refiner;
      //! If we set a FineSlowness model we will do the refinement otherwise we save some time by skipping it
      bool wantrefinement;
      //! The observed travel times for each source receiver pair
      jiba::rvec ObservedData;
      //! The calculator object
      TomographyCalculator Calculator;
      //!The implementation of the function to calculate the difference between observed and synthetic data
      virtual void
      ImplDataDifference(const jiba::rvec &Model, jiba::rvec &Diff);
      //! The implementation of the gradient calculation
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
      //! Set the model that contains the grid geometry and measurement and shot positions for the forward calculations
      void SetFineModelGeometry(const jiba::ThreeDSeismicModel &Model)
        {
          FineSlownessModel = Model;
          //set the coordinates of the refiner object according to the fine model
          //we want to use for the forward calculation
          Refiner.SetXCoordinates(FineSlownessModel.GetXCoordinates());
          Refiner.SetYCoordinates(FineSlownessModel.GetYCoordinates());
          Refiner.SetZCoordinates(FineSlownessModel.GetZCoordinates());
          //set the flag that we call the refiner in the forward and gradient calculation
          wantrefinement = true;
        }
      //! Set the geometry of the inversion model, the forward model will be determined by applying model refinement with FineSlownessModel
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
