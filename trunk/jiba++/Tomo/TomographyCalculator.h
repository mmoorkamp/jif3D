//============================================================================
// Name        : TomographyCalculator.h
// Author      : Apr 14, 2009
// Version     :
// Copyright   : 2009, mmoorkamp
//============================================================================


#ifndef TOMOGRAPHYCALCULATOR_H_
#define TOMOGRAPHYCALCULATOR_H_

#include "ThreeDSeismicModel.h"
#include "../Global/VecMat.h"
#include "modeling_seismic.h"

namespace jiba
  {
    /** \addtogroup Tomo Classes for seismic tomography */
    /* @{ */
    //! This class provides a convenient interface to the Podvin and Lecomte code and Bjoern's additions
    class TomographyCalculator
      {
    private:
      size_t nairlayers;
      jiba::GEOMETRY geo;
      jiba::GRID_STRUCT grid;
      jiba::DATA_STRUCT data;
      std::vector<jiba::RP_STRUCT> raypath;
      // Perform the dynamic allocation for the c-structures above
      void Allocate(const size_t ngrid, const size_t ndata, const size_t npos);
    public:
      TomographyCalculator();
      virtual ~TomographyCalculator();
      //! Return the raypath structure for the last forward modelling
      const std::vector<jiba::RP_STRUCT> &GetRayPath() const {return raypath;}
      //! Calculate the travel times for the given model
      /*! This is the core call to the forward modeling code. Has to be performed before LQDerivative can be called
       * @param Model The object containing the slowness distribution and measurement setup
       * @return The calculated travel times
       */
      rvec Calculate(const ThreeDSeismicModel &Model);
      //! Calculate the least-square derivative for the given model and data difference
      /*! For inversion we need the derivative of a least-squares objective function
       * with respect to the model parameters.
       * @param Model The object containing the slowness distribution and measurement setup, should be the same as the previous call to Calculate
       * @param Misfit The difference between observed and calculated data
       * @return The partial derivatives with respect to each model parameter
       */
      rvec LQDerivative(const ThreeDSeismicModel &Model, const rvec &Misfit);
      };
  /* @} */
  }

#endif /* TOMOGRAPHYCALCULATOR_H_ */
