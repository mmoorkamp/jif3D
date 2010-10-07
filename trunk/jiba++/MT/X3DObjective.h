//============================================================================
// Name        : X3DObjective.h
// Author      : Jul 10, 2009
// Version     : 
// Copyright   : 2009, mmoorkamp
//============================================================================


#ifndef X3DOBJECTIVE_H_
#define X3DOBJECTIVE_H_

#include "../ModelBase/ModelRefiner.h"
#include "../Global/FatalException.h"
#include "../Inversion/ObjectiveFunction.h"
#include "X3DModel.h"
#include "X3DMTCalculator.h"

namespace jiba
  {
    /** \addtogroup mtmodelling Forward modelling of magnetotelluric data */
    /* @{ */
    //! An objective function for 3D magnetotelluric data using x3D by Avdeev et al. as an engine
    class X3DObjective: public ObjectiveFunction
      {
    private:
      //! The object that calculates the synthetic data using x3d
      jiba::X3DMTCalculator Calculator;
      //! The skeleton for the conductivity model that contains geometry information etc.
      /*! The objective function simply copies the vector of values it receives for the current model
       * into this object, to obtain a well formed model for the synthetic forward and gradient calculation.
       * We therefore have to make sure that the total number of elements in the grid
       * matches the size of the model vector. Currently we cannot change the values for the background,
       * these are use as specified in ConductivityModel.
       *
       * There are a few restrictions on the conductivity model in order to ensure correct calculation of gradients. x3d
       * performs some "optimizations" that interfere with our approach for gradient calculation. Two adjacent layers
       * in the grid cannot have absolutely identical values, as x3d will then treat them as a single layer. Therefore at
       * least one grid cell must have a slightly different value from the layer above and below. Similarly a homogeneous
       * layer (with all grid cell values identical) cannot match the value for the background layer at this depth. In both
       * cases x3d will not yield field values for these areas in the .ema files that we need to calculate the gradient.
       * As the difference does not have to be very large, it is probably best to always have one grid cell in each layer
       * that differs from everything else by 0.1% or so.
       */
      X3DModel CoarseConductivityModel;
      //! Similarly to the seismic tomography case \see TomographyObjective we can specify a finer model for forward calculation
      X3DModel FineConductivityModel;
      //! Did we actually set a fine model for the forward calculation
      bool wantrefinement;
      //! The object performing the model refinement
      jiba::ModelRefiner Refiner;
      //! The vector of observed impedances real and imaginary part of the impedance elements for all stations at one frequency, then for the next frequency etc.
      jiba::rvec ObservedData;
      //! Calculate the difference between observed and synthetic data for a given conductivity model
      virtual void
      ImplDataDifference(const jiba::rvec &Model, jiba::rvec &Diff);
      //! The implementation of the gradient calculation
      virtual jiba::rvec ImplGradient(const jiba::rvec &Model,
          const jiba::rvec &Diff);
      //! So far transformations have no effect for MT data !
      virtual void SetDataTransformAction()
        {
          //TODO Implement transformation if necessary
        }
    public:
      //! Set the observed impedance data
      /*! We have to set the observed impedances for all sites and frequencies. They are stored
       * real and imaginary part separately, then Zxx, Zxy, Zux, Zyy, then by sites, then by frequency.
       * @param Data The observed data in the format described above.
       */
      void SetObservedData(const jiba::rvec &Data)
        {
          if (Data.empty())
            throw jiba::FatalException(
                "Cannot have empty observations in MT objective function.");
          ObservedData = Data;
        }
      //! Set a skeleton for the conductivity model that contains all information about cell sizes, site coordinates etc.
      /*! During the inversion we only copy conductivity values, so we have to store the
       * geometry and background information so we can form a complete model for the
       * forward calculation. The coarse model has to have the same number of cells as the model vector
       * that is passed to the gradient and forward calculations. It is used to establish the geometry
       * of the inversion model. For forward modeling this can be refined to the geometry of the FineConductivityModel
       * @param Model The skeleton object that contains all information about geometry, background and the location of the sites
       */
      void SetCoarseModelGeometry(const jiba::X3DModel &Model)
        {
          if (Model.GetConductivities().num_elements() == 0
              || Model.GetFrequencies().size() == 0)
            throw jiba::FatalException(
                "Cannot have empty frequencies or model in MT objective function.");
          CoarseConductivityModel = Model;
        }
      //! Set the geometry of a finely discretized model to ensure numerical precision
      /*! The geometry we want to use for the inversion might not be fine enough
       * to ensure a proper numerical calculation. In these cases we can pass the geometry
       * of a refined model to the objective function. Note that the resulting model
       * will have the grid interfaces from both the coarse model and the fine model
       * as described for the ModelRefiner class. The resulting model must still have
       * equal grid sizes in each of the horizontal directions. This is not checked here.
       * @param Model An object containing the refined grid coordinates, conductivity values, site locations
       *        and back are ignored.
       */
      void SetFineModelGeometry(const jiba::X3DModel &Model)
        {
          if (Model.GetConductivities().num_elements() == 0)
            throw jiba::FatalException(
                "Cannot have empty model in MT objective function.");
          FineConductivityModel = Model;
          //set the coordinates of the refiner object according to the fine model
          //we want to use for the forward calculation
          Refiner.SetXCoordinates(FineConductivityModel.GetXCoordinates());
          Refiner.SetYCoordinates(FineConductivityModel.GetYCoordinates());
          Refiner.SetZCoordinates(FineConductivityModel.GetZCoordinates());
          //set the flag that we call the refiner in the forward and gradient calculation
          wantrefinement = true;
        }
      X3DObjective();
      virtual ~X3DObjective();
      };
  /* @} */
  }

#endif /* X3DOBJECTIVE_H_ */
