//============================================================================
// Name        : Noise.h
// Author      : Oct 20, 2009
// Version     : 
// Copyright   : 2009, mmoorkamp
//============================================================================


#ifndef NOISE_H_
#define NOISE_H_

#include "VecMat.h"
#include <ctime>
#include <boost/random/lagged_fibonacci.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/bind.hpp>

namespace jiba
  {
    void AddNoise(jiba::rvec &Data, const double relerror,
        const double abserror)
      {

        boost::lagged_fibonacci607 generator(
            static_cast<unsigned int> (std::time(0)));
        const size_t ndata = Data.size();
        for (size_t i = 0; i < ndata; ++i)
          {
            const double error = std::max(std::abs(Data(i)) * relerror,
                abserror);
            boost::normal_distribution<> noise_dist(Data(i), error);
            boost::variate_generator<boost::lagged_fibonacci607&,
                boost::normal_distribution<> > noise(generator, noise_dist);
            Data(i) = noise();
          }

      }
  }
#endif /* NOISE_H_ */
