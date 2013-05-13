//============================================================================
// Name        : OneDMTCalculator.h
// Author      : 13 Jun 2012
// Version     : 
// Copyright   : 2012, mm489
//============================================================================

#ifndef ONEDMTCALCULATOR_H_
#define ONEDMTCALCULATOR_H_

#include "../Global/VecMat.h"
#include "X3DModel.h"

namespace jiba
  {
    /** \addtogroup mtmodelling Forward modelling of magnetotelluric data */
    /* @{ */
    //! This class uses the Cagniard algorithm to calculate the MT response of a 1D layered halfspace
    /*! This is an implementation of the Cagniard algorithm
     * to calculate magnetotelluric impedances for a layered half-space. It also
     * implements the calculation of the derivative of a least-squares objective
     * function using an adjoint approach and the equations in Avdeeva, 2006.
     */
    class OneDMTCalculator
      {
    public:
      //! This type definition is necessary so that ThreeDModelObjective can correctly deduce the native type for a model object for this class
      typedef X3DModel ModelType;
    private:
      jiba::cmat alpha;
      jiba::cmat gammakj;
      jiba::cmat gammaj;
      jiba::cvec Z;
    public:
      //! Given a model of layer thicknesses and conductivities we calculate MT impedances
      rvec Calculate(const ModelType &Model);
      //! Given a model and a misfit vector, calculate the derivative of a least-squares objective function using an adjoint approach
      rvec LQDerivative(const ModelType &Model, const rvec &Misfit);
      OneDMTCalculator();
      virtual ~OneDMTCalculator();
      };
  /* @} */
  } /* namespace jiba */
#endif /* ONEDMTCALCULATOR_H_ */