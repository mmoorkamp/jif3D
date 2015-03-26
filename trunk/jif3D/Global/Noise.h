//============================================================================
// Name        : Noise.h
// Author      : Oct 20, 2009
// Version     : 
// Copyright   : 2009, mmoorkamp
//============================================================================

#ifndef NOISE_H_
#define NOISE_H_

#include "FatalException.h"
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

namespace jif3D
  {
    /** \addtogroup util General utility routines */
    /* @{ */

    //! Add random gaussian noise to a vector of data
    /*! When inverting synthetic data, we usually want to add
     * noise to it to simulate some of the properties of real data. This function
     * takes a vector of real data values and adds noise to it. The variance
     * can be specified both as a relative or absolute value. If both are specified the maximum is taken.
     * @param Data The vector of data values, contains the noisy data afterwards
     * @param relerror The relative error for each datum, e.g. 0.02 corresponds to 2%
     * @param abserror A vector of minimum absolute errors for each datum in the same units as the data vector, has to have the same length as the data vector
     */
    inline void AddNoise(jif3D::rvec &Data, const double relerror, jif3D::rvec &abserror)
      {
        //create a random number generator object without specifying the distribution
        boost::lagged_fibonacci607 generator(static_cast<unsigned int>(std::time(0)));
        const size_t ndata = Data.size();
        //check that we have an error for each datum
        if (Data.size() != abserror.size())
          throw jif3D::FatalException(
              "Data and Error vectors do not have the same size !");
        //create an error estimate for each datum
        for (size_t i = 0; i < ndata; ++i)
          {
            const double error = std::max(std::abs(Data(i)) * relerror, abserror(i));
            //construct a Gaussian distribution based on mean and standard deviation
            boost::normal_distribution<> noise_dist(Data(i), error);
            boost::variate_generator<boost::lagged_fibonacci607&,
                boost::normal_distribution<> > noise(generator, noise_dist);
            Data(i) = noise();
          }

      }

    //! Add random gaussian noise to a vector of data
    /*! When inverting synthetic data, we usually want to add
     * noise to it to simulate some of the properties of real data. This function
     * takes a vector of real data values and adds noise to it. The variance
     * can be specified both as a relative or absolute value. If both are specified the maximum is taken.
     * @param Data The vector of data values, contains the noisy data afterwards
     * @param relerror The relative error for each datum, e.g. 0.02 corresponds to 2%
     * @param abserror The minimum absolute error for all data in the same units as the data vector, this value is applied globally
     */
    inline void AddNoise(jif3D::rvec &Data, const double relerror, const double abserror =
        0)
      {
        jif3D::rvec AbsErrorVec(Data.size(), abserror);
        AddNoise(Data, relerror, AbsErrorVec);
      }

    //! This function provides a simple, convenient way to assign an error to data for inversion based on actual errors and error floor
    /*! When inverting data we often have to use a minimum error based on some noise floor estimate.
     * Also when inverting noise free synthetic data, we still have
     * to assign some data covariance for inversion, usually we use a relative error.
     * This function provides such functionality and
     * takes care of issues with zero or very small data.
     * @param Data The vector of data elements
     * @param relerror The minimum relative error for each datum
     * @param absmin The absolute minimum data value considered for error calculation, this reduced the influence of very small data
     * @return The vector of error estimates
     */
    inline jif3D::rvec ConstructError(const jif3D::rvec &Data,
        const jif3D::rvec &DataError, const double relerror, const double absmin = 0.0)
      {
        //check for reasonable relative error value
        if (relerror <= 0.0)
          {
            throw jif3D::FatalException("Specifiying relative error <= 0 is not valid!");
          }
        //check for reasonable absolute error value
        if (absmin <= 0.0)
          {
            throw jif3D::FatalException("Specifiying absolute error <= 0 is not valid!");
          }
        const size_t ndata = Data.size();
        //create objects for the misfit and a very basic error estimate
        jif3D::rvec Error(ndata, 0.0);
        for (size_t i = 0; i < ndata; ++i)
          {
            double minerr = std::max(std::abs(Data(i)) * relerror, absmin);
            Error(i) = std::max(DataError(i), minerr);
            assert(Error(i) > 0.0);
          }
        return Error;
      }
    //! Assign errors to noise free synthetic MT data for inversion purposes
    /*! When inverting synthetic data for which we do not have any noise information, we still have
     * to assign some data covariance for inversion, usually we use a relative error. This function, in contrast
     * to ConstructError, is particularly geared towards MT data. We use a percentage of the Berdichevskyi invariant
     * for a given frequency and site as the error estimate for all other tensor elements at that
     * frequency and site. This avoids problems with small diagonal elements.
     * @param Data A vector containing the MT data, we assume that 8 consecutive real numbers form one complex impedance tensor
     * @param relerror The relative error of the maximum tensor element
     * @return The vector of error estimates
     */
    inline jif3D::rvec ConstructMTError(const jif3D::rvec &Data, const double relerror)
      {
        if (relerror <= 0.0)
          {
            throw jif3D::FatalException("Specifiying relative error <= 0 is not valid!");
          }
        const size_t ndata = Data.size();
        jif3D::rvec DataError(ndata, 0.0);
        const size_t ntensorelem = 8;
        if ((Data.size() % ntensorelem) != 0)
          {
            throw jif3D::FatalException(
                "MT Data vector size is not an integer multiple of 8!");
          }
        const size_t ntensor = ndata / ntensorelem;
        for (size_t i = 0; i < ntensor; ++i)
          {
            //compute the real and imaginary parts of the berdichevskyi invariant
            //Berd = 0.5 * (Zxy - Zyx)
            double berdreal = 0.5
                * (Data(i * ntensorelem + 2) - Data(i * ntensorelem + 4));
            double berdimag = 0.5
                * (Data(i * ntensorelem + 3) - Data(i * ntensorelem + 5));
            // we assume the absolute value of the invariant as a reference threshold
            //for the error calculation
            double berdabs = sqrt(berdreal * berdreal + berdimag * berdimag);
            std::fill_n(DataError.begin() + i * ntensorelem, ntensorelem,
                berdabs * relerror);
          }
        return DataError;
      }
  /* @} */
  }
#endif /* NOISE_H_ */
