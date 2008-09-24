//============================================================================
// Name        : DepthWeigthing.cpp
// Author      : Sep 18, 2008
// Version     :
// Copyright   : 2008, mmoorkamp
//============================================================================

#include <fstream>
#include <iostream>
#include "DepthWeighting.h"
#include "../Inversion/LinearInversion.h"

namespace jiba
  {

    double WeightingTerm::operator()(const double z, const double z0) const
      {
        return pow(z + z0, n);
      }
    double WeightingTerm::deriv(const double z, const double z0) const
      {
        return n * pow(z + z0, n-1);
      }
    double WeightingTerm::average(const double z1, const double z2, const double z0) const
      {
        return 1./((z2-z1)*(n+1)) * (pow(z2 + z0, n+1) - pow(z1 + z0, n+1));
      }
    /*!
     * @param SensProfile A vector with nz elements that contains the sensitivity (for a single measurement) with depth
     * @param ZSizes The size of each cell in z-direction
     * @param WeightFunction The function object that gives the form of the weighting terms and it derivatives
     * @return An estimate of z0 for which the weighting function matches the sensitivity decay most closely
     */
    double FitZ0(const jiba::rvec &SensProfile,
        const ThreeDModelBase::t3DModelDim &ZSizes, const jiba::WeightingTerm &WeightFunction)
      {


        const size_t zsize = SensProfile.size();
        jiba::rvec profile(SensProfile), zvalues(zsize);
        //we need the depth to the interfaces
        partial_sum(ZSizes.begin(), ZSizes.end(),
            zvalues.begin());
        //the sensitivities can be negative, so we examine the absolute value
        transform(profile.begin(),profile.end(),profile.begin(),std::abs<double>);


        //then we normalize so the largest element (usually the top cell) has value 1
        transform(profile.begin(), profile.end(), profile.begin(),
            boost::bind(std::divides<double>(), _1, *std::max_element(
                    profile.begin(), profile.end())));
        // we guess a value for z0, we take it slightly larger than the bottom of the model
        double startz = zvalues(zsize - 1) * 1.1;

        //setup the "inversion"
        const size_t ndata = zvalues.size();
        jiba::rmat sens(zsize, 1);
        jiba::rvec weights(1);
        jiba::rvec InvModel(1), DeltaModel(1);
        weights(0) = 1.0;
        InvModel(0) = startz;
        const size_t iterations = 100;
        jiba::rvec error(ndata), delta(ndata), calculated(ndata);
        std::fill_n(error.begin(), zsize, 1.0);
        const double evalthresh = 1e-6;
        std::ofstream outfile("match.out");
        double stepsize = 1e6;
        size_t i = 0;
        //while we are still making progress, but don't use too many iterations
        while (stepsize> 0.1 && i < iterations)
          {
            //get the predicted data
            for (size_t j = 0; j < zsize; ++j)
              {
                calculated(j) = WeightFunction(zvalues(j),InvModel(0));
              }
            double maximum = *std::max_element(calculated.begin(),
                calculated.end());
            //normalize the predicted data just like the sensitivities
            transform(calculated.begin(), calculated.end(), calculated.begin(),
                boost::bind(std::divides<double>(), _1, maximum));
            //calculate the derivatives and the misfit
            for (size_t j = 0; j < zsize; ++j)
              {
                sens(j, 0) = WeightFunction.deriv(zvalues(j), InvModel(0)) / maximum;
                delta(j) = profile(j) - calculated(j);
              }
            for (size_t j = 0; j < zsize; ++j)
              {
                outfile << zvalues(j) << " " << profile(j) << " "
                << calculated(j) << std::endl;
              }
            // do one step of the inversion
            jiba::DataSpaceInversion()(sens, delta, weights, error, evalthresh,
                0.0, DeltaModel);
            //and adjust the model (no fancy line search here)
            stepsize = boost::numeric::ublas::norm_2(DeltaModel);
            InvModel -= DeltaModel;

            outfile << std::endl << std::endl;
            ++i;
          }

        return InvModel(0);
      }

    /*! Construct a weighting vector that counteracts the decay of the sensitivities of gravity data
     *  to facilitate inversion.
     * @param XSizes The sizes of the model cells in x-direction in m
     * @param YSizes The sizes of the model cells in y-direction in m
     * @param ZSizes The sizes of the model cells in z-direction in m
     * @param z0 The scaling depth
     * @param WeightVector The resulting weight vector
     *
     * The resulting weight vector has only nz elements, where nz is the number of vertical layers
     * the inversion algorithm has to take care to apply these values to the appropriate matrix elements
     */
    void ConstructDepthWeighting(const ThreeDModelBase::t3DModelDim &XSizes,
        const ThreeDModelBase::t3DModelDim &YSizes,
        const ThreeDModelBase::t3DModelDim &ZSizes, const double z0,
        rvec &WeightVector, const jiba::WeightingTerm &WeightFunction)
      {
        const size_t nz = ZSizes.size();

        if (WeightVector.size() != nz)
          {
            WeightVector.resize(nz);
          }
        double currbottom = 0.0;
        for (size_t i = 0; i < nz; ++i)
          {
            currbottom += ZSizes[i];
            WeightVector(i) = WeightFunction(currbottom,z0);

          }
        //normalize
        std::transform(WeightVector.begin(), WeightVector.end(),
            WeightVector.begin(), boost::bind(std::divides<double>(), _1,
                *std::max_element(WeightVector.begin(), WeightVector.end())));

      }

  }
