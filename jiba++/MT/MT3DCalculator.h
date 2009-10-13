//============================================================================
// Name        : MT3DCalculator.h
// Author      : Jul 6, 2009
// Version     : 
// Copyright   : 2009, mmoorkamp
//============================================================================


#ifndef MT3DCALCULATOR_H_
#define MT3DCALCULATOR_H_

#include "../Global/VecMat.h"
#include "ThreeDMTModel.h"

namespace jiba
  {
    /** \addtogroup mtmodelling Forward modelling of magnetotelluric data */
    /* @{ */
    //! This is the base class for the implementation of the 3D MT forward calculation
    /*! For our joint inversion we use x3D by Avdeev et al. to calculate impedances
     * for MT. However, this code is not generally available so we provide this abstract
     * base class that defines the interface for the MT forward calculation. By deriving
     * from this class and implementing the two virtual functions CalculateImpl and
     * LQDerivativeImpl any other implementation can be plugged into the joint inversion
     * with only minor changes to the code.
     */
    class MT3DCalculator
      {
    private:
      //! Given a conductivity model, calculate a vector of impedances
      virtual rvec CalculateImpl(const ThreeDMTModel &Model) = 0;
      //! Given a conductivity model and the misfit for each datum, calculate the derivative of the objective function with respect to the model parameters.
      virtual rvec LQDerivativeImpl(const ThreeDMTModel &Model,
          const rvec &Misfit) = 0;
    public:
      //! Calculate synthetic impedances for a 3D MT model, forwards call to CalculateImpl
      rvec Calculate(const ThreeDMTModel &Model);
      //! Given a conductivity model and the misfit for each datum, calculate the derivative of the objective function with respect to the model parameters.
      rvec LQDerivative(const ThreeDMTModel &Model, const rvec &Misfit);
      MT3DCalculator();
      virtual ~MT3DCalculator();
      };
  /* @} */
  }

#endif /* MT3DCALCULATOR_H_ */
