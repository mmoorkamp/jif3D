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
