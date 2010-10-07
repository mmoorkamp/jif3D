//============================================================================
// Name        : GravityObjective.h
// Author      : Apr 16, 2009
// Version     :
// Copyright   : 2009, mmoorkamp
//============================================================================


#ifndef GRAVITYOBJECTIVE_H_
#define GRAVITYOBJECTIVE_H_

#include "../Inversion/ObjectiveFunction.h"
#include "ThreeDGravityCalculator.h"
#include <boost/shared_ptr.hpp>

namespace jiba
  {
    /** \addtogroup gravity Gravity forward modeling, display and inversion */
    /* @{ */

    //! This class defines the objective function for both scalar and FTG gravity data
    /*! We have a single objective function class for scalar and tensor gravity data, which
     * type the data vector contains is selected in the constructor and scalar data by default.
     *
     * We can choose to minimize derived quantities like an invariant for FTG through the data transformation
     * mechanism that is described in the base class ObjectiveFunction.
     *
     * Before we can use an object of this type we have to set the model geometry correctly. This can
     * be achieved through the SetModelGeometry function. The Model passed to this function has to
     * have an appropriate size for the model vector that the objective function gets in the CalcMisfit call.
     * If the size of the model vector matches the size of the gridded part it is interpreted to represent the
     * grid. If matches the size of the gridded part plus the background densities, both the grid and the background
     * are updated. Also the measurement positions in the model have to match the number of observed data.
     */
    class GravityObjective: public ObjectiveFunction
      {
    private:
      //! The pointer to the forward calculation object, assigen by the constructor
      boost::shared_ptr<jiba::ThreeDGravityCalculator> Calculator;
      //! The observed gravity data, its type (scalar/FTG) depends on the constructor
      jiba::rvec ObservedData;
      //! The object containing the geometry information of the model, the densities are overwritten at each misfit calculation
      jiba::ThreeDGravityModel DensityModel;
      //Calculate the difference between observed and synthetic data from the current model
      virtual void ImplDataDifference(const jiba::rvec &Model, jiba::rvec &Diff);
      // Calculate the model gradient
      virtual jiba::rvec ImplGradient(const jiba::rvec &Model, const  jiba::rvec &Diff);
      //When we set a transform in the base class, we have to pass this transform to the forward calculation class
      virtual void SetDataTransformAction()
        {
          Calculator->SetDataTransform(GetDataTransform());
        }
    public:
      //! Set the observed data
      /*! This function takes a simple vector of observations that can represent either scalar or FTG data.
       * How this vector is interpreted depends with which parameters the objective function object was created.
       * If it represents FTG data we expect 9 elements per measurement site.
       * @param Data
       */
      void SetObservedData(const jiba::rvec &Data)
        {
          ObservedData = Data;
        }
      //! Set the model object that contains the geometry information, cell sizes, background and measurement positions
      void SetModelGeometry(const jiba::ThreeDGravityModel &Model)
        {
          DensityModel = Model;
        }
      //! During construction we have to choose the data type and whether we want to use CUDA for forward calculations
      /*! The constructor takes two arguments
       * @param ftg Do we minimize FTG data, if yes we expect 9 elements per measurement, default false
       * @param cuda Use CUDA for forward calculations, default false
       */
      GravityObjective(bool ftg = false, bool cuda = false);
      virtual ~GravityObjective();
      };
  /* @} */
  }

#endif /* GRAVITYOBJECTIVE_H_ */
