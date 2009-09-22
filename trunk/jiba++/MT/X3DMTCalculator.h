//============================================================================
// Name        : X3DMTCalculator.h
// Author      : Jul 7, 2009
// Version     : 
// Copyright   : 2009, mmoorkamp
//============================================================================


#ifndef X3DMTCALCULATOR_H_
#define X3DMTCALCULATOR_H_

#include "MT3DCalculator.h"
#include "X3DModel.h"

namespace jiba
  {
    /** \addtogroup mtmodelling Forward modelling of magnetotelluric data */
    /* @{ */
    class X3DMTCalculator
      {
    private:
      //! Create a unique name for each object, calculation type and frequency so that we can write to different directories and execute in parallel
      std::string MakeUniqueName(X3DModel::ProblemType Type,
          const size_t FreqIndex);
    public:
      //! Given a conductivity model, calculate a vector of impedances
      rvec Calculate(const X3DModel &Model);
      //! Given a conductivity model and the misfit for each datum, calculate the derivative of the objective function with respect to the model parameters.
      rvec LQDerivative(const X3DModel &Model, const rvec &Misfit);

      X3DMTCalculator();
      virtual ~X3DMTCalculator();
      };
  /* @} */
  }

#endif /* X3DMTCALCULATOR_H_ */
