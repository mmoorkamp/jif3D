//============================================================================
// Name        : ModEMCalculator.h
// Author      : 14 Jan 2016
// Version     : 
// Copyright   : 2016, mm489
//============================================================================

#ifndef MT_MODEMCALCULATOR_H_
#define MT_MODEMCALCULATOR_H_

#include "../Global/Serialization.h"
#include "../Global/VecMat.h"
#include "ThreeDMTModel.h"

#include <limits>

namespace jif3D
  {

    class ModEMCalculator
      {
    public:
      typedef ThreeDMTModel ModelType;
      ModEMCalculator();
      virtual ~ModEMCalculator();
      //! Given a conductivity model, calculate a vector of impedances
      /*! For a conductivity model given by the input parameter Model, we calculate the synthetic magnetotelluric data. When compiled with
       * an appropriate compiler the calculation is run in parallel for each frequency. We return the synthetic data as a real valued vector.
       * The ordering is \f$Re(Z_xx),Im(Z_xx),Re(Z_xy),\ldots,Im(Z_yy)\f$ for the first frequency for all sites, then second frequency for all sites etc.
       *
       * @param Model The description of the conductivity model including sites locations and frequencies.
       * @param minfreqindex The index of the first frequency for which to calculate the gradient
       * @param maxfreqindex The index one larger than the index of the last frequency for which to calculate the gradient (C++ loop convention)
       * @return The synthetic MT data in the format described above.
       */
      rvec Calculate(const ModelType &Model, size_t minfreqindex = 0,
          size_t maxfreqindex = std::numeric_limits<size_t>::max());
      //! Given a conductivity model and the misfit for each datum, calculate the derivative of the objective function with respect to the model parameters.
      /*! We use an adjoint approach to calculate the gradient of the objective functions with respect to the model parameters. As this approach requires
       * some of the fields from the forward calculation, the gradient will only be correct if the function Calculate of the same object has been called for
       * the same model beforehand. It is safe to calculate different models with separate objects between those calls.
       * @param Model The description of the conductivity model. Has to be the same as for the previous call to calculate.
       * @param Misfit The data misfit associated with the model.
       * @param minfreqindex The index of the first frequency for which to calculate the gradient
       * @param maxfreqindex The index one larger than the index of the last frequency for which to calculate the gradient (C++ loop convention)
       * @return The gradient of the objective function with respect to the model parameters for the given model. The storage ordering is identical to X3DModel.
       */
      rvec LQDerivative(const ModelType &Model, const rvec &Misfit, size_t minfreqindex =
          0, size_t maxfreqindex = std::numeric_limits<size_t>::max());
      };

  } /* namespace jif3D */

#endif /* MT_MODEMCALCULATOR_H_ */
