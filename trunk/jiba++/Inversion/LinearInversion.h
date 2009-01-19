//============================================================================
// Name        : LinearInversion.h
// Author      : Max Moorkamp
// Version     :
// Copyright   : 2008, MM
//============================================================================


#ifndef LINEARINVERSION_H_
#define LINEARINVERSION_H_

#include "../Global/VecMat.h"

/*! \file LinearInversion.h
 * Linear inversion routines for general use.
 */
namespace jiba
  {
    /** \addtogroup inversion General routines for inversion
     * These are all general inversion classes that implement standard
     * inversion methods. These classes know nothing about the type of
     * data or specific method involved, but purely implement the mathematical
     * operations. It is up to the inversion program or more specific classes
     * to provide a link and accommodate the peculiarities of the geophysical
     * measurements we want to invert.
     */
    /* @{ */

    //! Invert data described by a linear relationship using the dataspace method
    class DataSpaceInversion
      {
    public:
      void operator()(rmat &Sensitivities, rvec &Data, const rvec &WeightVector,
          const rvec &DataError, const double lambda, rvec &InvModel);
      };

    //! Invert data described by a linear relationship using the model space method
    class ModelSpaceInversion
      {
    public:
      void operator()(rmat &Sensitivities, rvec &Data, const rvec &WeightVector,
          const rvec &DataError, const double lambda,
          rvec &InvModel);
      };
    //! Calculate a model update using Quasi-Newton Inversion
    class QuasiNewtonInversion
      {
    public:
      void operator()(rmat &Sensitivities, rvec &Data, const rvec &WeightVector,
          const rvec &DataError, const double lambda, rvec &InvModel);
      };
  /* @} */
  }

#endif /* LINEARINVERSION_H_ */
