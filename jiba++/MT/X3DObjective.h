//============================================================================
// Name        : X3DObjective.h
// Author      : Jul 10, 2009
// Version     : 
// Copyright   : 2009, mmoorkamp
//============================================================================


#ifndef X3DOBJECTIVE_H_
#define X3DOBJECTIVE_H_

#include "../Inversion/ObjectiveFunction.h"
#include "X3DModel.h"
#include "X3DMTCalculator.h"

namespace jiba
  {
    /** \addtogroup mtmodelling Forward modelling of magnetotelluric data */
    /* @{ */
    class X3DObjective: public ObjectiveFunction
      {
    private:
      //The object that calculates the synthetic data using x3d
      jiba::X3DMTCalculator Calculator;
      //The skeleton for the conductivity model that contains geometry information etc.
      //the objective function simply copies the vector of value it receives for the current model
      //into this object, to obtain a well formed model for the synthetic forward and gradient calculation
      X3DModel ConductivityModel;
      //the vector of observed impedances real and imaginary part of the impedance elements
      //for all stations at one frequency, then for the next frequency etc.
      jiba::rvec ObservedData;
      //calculate the difference between observed and synthetic data for a given conductivity model
      virtual void
      ImplDataDifference(const jiba::rvec &Model, jiba::rvec &Diff);
      //The implementation of the gradient calculation
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
       * real and imaginary part seperately, then Zxx, Zxy, Zux, Zyy, then by sites, then by frequency.
       * @param Data The observed data in the format described above.
       */
      void SetObservedData(const jiba::rvec &Data)
        {
          ObservedData = Data;
        }
      //! Set a skeleton for the conductivity model that contains all information about cell sizes, site coordinates etc.
      /*! During the inversion we only copy conductivity values, so we have to store the
       * geometry and background information so we can form a complete model for the
       * forward calculation.
       * @param Model The skeleton object that contains all information about geometry and background
       */
      void SetModelGeometry(const jiba::X3DModel &Model)
        {
          ConductivityModel = Model;
        }
      X3DObjective();
      virtual ~X3DObjective();
      };
  /* @} */
  }

#endif /* X3DOBJECTIVE_H_ */
