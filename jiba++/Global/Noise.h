//============================================================================
// Name        : Noise.h
// Author      : Oct 20, 2009
// Version     : 
// Copyright   : 2009, mmoorkamp
//============================================================================


#ifndef NOISE_H_
#define NOISE_H_

#include "VecMat.h"
#include "NumUtil.h"
#include <ctime>
#include <boost/random/lagged_fibonacci.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/bind.hpp>

/*! \file Noise.h
 * General routines for adding noise to observed data or constructing error estimates
 * for the inversion of noisy data.
 */

namespace jiba
  {
    //! Add random gaussian noise to a vector of data
    /*! When inverting synthetic data, we usually want to add
     * noise to it to simulate some of the properties of real data. This function
     * takes a vector of real data values and adds noise to it. The variance
     * can be specified both as a relative or absolute value. If both are specified the maximum is taken.
     * @param Data The vector of data values, contains the noisy data afterwards
     * @param relerror The relative error for each datum, e.g. 0.02 corresponds to 2%
     * @param abserror The minimum absolute error for each datum in the same units as the data vector
     */
    inline void AddNoise(jiba::rvec &Data, const double relerror,
        const double abserror = 0)
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

    //! This function provides a simple, convenient way to assign an error to synthetic data for inversion
    /*! When inverting synthetic data for which we do not have any noise information, we still have
     * to assign some data covariance for inversion, usually we use a relative error.
     * This function provides such functionality and
     * takes care of issues with zero or very small data.
     * @param Data The vector of data elements
     * @param errorlevel The minimum relative error for each datum
     * @param absmin The absolute minimum data value considered for error calculation, this reduced the influence of very small data
     * @return The vector of error estimates
     */
    inline jiba::rvec ConstructError(const jiba::rvec &Data, const double relerror,
        const double absmin = 0.0)
      {
        assert(relerror >= 0.0);
        assert(absmin >= 0.0);
        const size_t ndata = Data.size();
        //create objects for the misfit and a very basic error estimate
        jiba::rvec DataError(ndata);
        for (size_t i = 0; i < ndata; ++i)
          {
            DataError(i) = std::max(std::abs(Data(i)) * relerror, absmin);
            assert(DataError(i) > 0.0);
          }
        return DataError;
      }
    //! Assign errors to noise free synthetic MT data for inversion purposes
    /*! When inverting synthetic data for which we do not have any noise information, we still have
     * to assign some data covariance for inversion, usually we use a relative error. This function, in contrast
     * to ConstructError, is particularly geared towards MT data. We use a percentage of the maximum tensor element
     * for a given frequency and site as the error estimate for all other tensor elements at that
     * frequency and site. This avoids problems with small diagonal elements.
     * @param Data A vector containing the MT data, we assume that 8 consecutive real numbers form one complex impedance tensor
     * @param The relative error of the maximum tensor element
     * @return The vector of error estimates
     */
    inline jiba::rvec ConstructMTError(const jiba::rvec &Data, const double relerror)
      {
        assert(relerror >= 0.0);
        const size_t ndata = Data.size();
        jiba::rvec DataError(ndata);
        const size_t ntensorelem = 8;
        assert((Data.size() % ntensorelem) == 0);
        const size_t ntensor = ndata / 8;
        for (size_t i = 0; i < ntensor; ++i)
          {
            //find the maximum tensor element
            const double maxdata = std::abs(*std::max_element(Data.begin() + i
                * ntensorelem, Data.begin() + (i + 1) * ntensorelem,
                jiba::absLess<double, double>()));
            assert(maxdata > 0);
            std::fill_n(DataError.begin() + i * ntensorelem, ntensorelem,
                maxdata * relerror);
          }
        return DataError;
      }
  }
#endif /* NOISE_H_ */
