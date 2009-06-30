//============================================================================
// Name        : ConstructError.h
// Author      : May 22, 2009
// Version     :
// Copyright   : 2009, mmoorkamp
//============================================================================


#ifndef CONSTRUCTERROR_H_
#define CONSTRUCTERROR_H_

#include <cassert>
#include "../Global/VecMat.h"

namespace jiba
  {
    //! This function provides a simple, convenient way to assign an error to synthetic data for inversion
    /*! When inverting synthetic data for which we do not have any noise information, we still have
     * to assign some data covariance for inversion, usually we use a relative error.
     * This function provides such functionality and
     * takes care of issues with zero or very small data.
     * @param Data The vector of data elements
     * @param errorlevel The minimum relative error for each datum
     * @return The vector of error estimates
     */
    jiba::rvec ConstructError(const jiba::rvec &Data, const double errorlevel)
      {
        assert(errorlevel > 0.0);
        const size_t ndata = Data.size();
        const double maxdata = fabs(*std::max_element(Data.begin(), Data.end(),
            jiba::absLess<double, double>()));
        //create objects for the misfit and a very basic error estimate
        jiba::rvec DataError(ndata);
        for (size_t i = 0; i < ndata; ++i)
          {
            DataError(i) = std::max(std::fabs(Data(i)), 1e-2 * maxdata)
                * errorlevel;
            assert(DataError(i) > 0.0);
          }
        return DataError;
      }
  }
#endif /* CONSTRUCTERROR_H_ */
