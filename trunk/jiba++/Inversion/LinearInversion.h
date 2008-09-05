//============================================================================
// Name        : gravgrid.cpp
// Author      : Max Moorkamp
// Version     :
// Copyright   : 2008, MM
//============================================================================


#ifndef LINEARINVERSION_H_
#define LINEARINVERSION_H_

#include "../Global/VecMat.h"

namespace jiba
  {
    /** \addtogroup inversion General routines for inversion */
    /* @{ */

    //! Invert data described by a linear relationship using the dataspace method
    class DataSpaceInversion
      {
    private:
      jiba::rmat FilteredSens;
    public:
      const jiba::rmat &GetFilteredSens()
        {
          return FilteredSens;
        }
      void operator()(rmat &Sensitivities, rvec &Data, const rvec &WeightVector,
          const rvec &DataError, const double evalthresh, const double lambda, rvec &InvModel);
      };

    class ModelSpaceInversion
      {
    private:
      jiba::rmat FilteredSens;
    public:
      const jiba::rmat &GetFilteredSens()
        {
          return FilteredSens;
        }
      void operator()(rmat &Sensitivities, rvec &Data, const rvec &WeightVector,
          const rvec &DataError, const double evalthresh,const double lambda, rvec &InvModel);
      };
  /* @} */
  }

#endif /* LINEARINVERSION_H_ */
