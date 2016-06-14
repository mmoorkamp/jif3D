//============================================================================
// Name        : MT3DCalculator.cpp
// Author      : Jul 6, 2009
// Version     : 
// Copyright   : 2009, mmoorkamp
//============================================================================


#include "MT3DCalculator.h"

namespace jif3D
  {
    /*! Calculate the response of a synthetic MT model. All information
     * about frequencies, site locations etc. is contained in the Model
     * Object.
     * @param Model The conductivity model we want to calculate the impedances for
     * @return A real vector of impedance values with the ordering \f$Re(Z_xx),Im(Z_xx),Re(Z_xy),\ldots,Im(Z_yy)\f$ for the first frequency for all sites, then second frequency for all sites etc.
     */
    rvec MT3DCalculator::Calculate(const ThreeDMTModel &Model)
      {
        return CalculateImpl(Model);
      }

    /*! Calculate the derivative of a least-squares objective function of the impedance elements.
     * @param Model The conductivity model we want to calculate the derivative for
     * @param Misfit The misfit for each impedance element in the same storage ordering as the synthetic data returned by Calculate
     * @return The real gradient vector of the objective function with the same storage ordering as the model
     */
    rvec MT3DCalculator::LQDerivative(const ThreeDMTModel &Model,
        const rvec &Misfit)
      {
        return LQDerivativeImpl(Model, Misfit);
      }

    MT3DCalculator::MT3DCalculator()
      {

      }

    MT3DCalculator::~MT3DCalculator()
      {

      }

  }
