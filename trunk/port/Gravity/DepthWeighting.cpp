//============================================================================
// Name        : DepthWeigthing.cpp
// Author      : Sep 18, 2008
// Version     :
// Copyright   : 2008, mmoorkamp
//============================================================================

#include <cmath>
#include <fstream>
#include <iostream>
#include "DepthWeighting.h"
#include "../Inversion/LinearInversion.h"

namespace jif3D
  {
    /*! Calculate the value of the function f
     * @param z The depth value
     * @param z0 The parametric depth value
     * @return \f$ (z + z_0)^n \f$
     */
    double WeightingTerm::operator()(const double z, const double z0) const
      {
        return pow(z + z0, n);
      }
    /*! Calculate the derivative of the function f at  given z and z0
     * @param z The depth value
     * @param z0 The parametric depth value
     * @return The derivative of the depth weighting function
     */
    double WeightingTerm::deriv(const double z, const double z0) const
      {
        return n * pow(z + z0, n - 1);
      }
    /*! Calculate the average value of f in the interval [z1,z2]
     * @param z1 The lower bound of the interval
     * @param z2 The upper bound of the interval
     * @param z0 The parametric depth value
     * @return The average value of f
     */
    double WeightingTerm::average(const double z1, const double z2, const double z0) const
      {
        return 1. / ((z2 - z1) * (n + 1)) * (pow(z2 + z0, n + 1) - pow(z1 + z0, n + 1));
      }
    /*!
     * @param SensProfile A vector with nz elements that contains the sensitivities (for a single measurement) for each cell with depth
     * @param ZSizes The size of each cell in z-direction
     * @param WeightFunction The function object that gives the form of the weighting terms and it derivatives, e.g. an instance of the class WeightingTerm
     * @return An estimate of z0 for which the weighting function matches the sensitivity decay most closely
     */
    double FitZ0(const jif3D::rvec &SensProfile,
        const ThreeDModelBase::t3DModelDim &ZSizes,
        const jif3D::WeightingTerm &WeightFunction)
      {

        const size_t zsize = SensProfile.size();
        jif3D::rvec profile(SensProfile), zvalues(zsize);
        //we need the depth to the interfaces
        partial_sum(ZSizes.begin(), ZSizes.end(), zvalues.begin());
        //the sensitivities can be negative, so we examine the absolute value
        for ( double & value : profile )
          {
            value = std::abs(value);
          }
        //then we normalize so the largest element (usually the top cell) has value 1
        const double maxel = *std::max_element(profile.begin(), profile.end());
        transform(profile.begin(), profile.end(), profile.begin(), [maxel] (double val)
          { return val/maxel;});
        // we guess a value for z0, we take it slightly larger than the bottom of the model
        //it shouldn't be on a cell boundary to avoid numerical problems
        double startz = zvalues(zsize - 1) * 1.1;

        //setup the "inversion"
        const size_t ndata = zvalues.size();
        jif3D::rmat sens(zsize, 1);
        jif3D::rvec weights(1);
        jif3D::rvec InvModel(1), DeltaModel(1);
        weights(0) = 1.0;
        InvModel(0) = startz;
        DeltaModel(0) = 0.0;
        const size_t iterations = 100;
        jif3D::rvec error(ndata), delta(ndata), calculated(ndata);
        std::fill_n(error.begin(), zsize, 1.0);

        std::ofstream outfile("match.out");
        double stepsize = 1e6;
        size_t i = 0;
        //while we are still making progress, but don't use too many iterations
        while (stepsize > 0.01 && i < iterations)
          {
            //get the predicted data
            for (size_t j = 0; j < zsize; ++j)
              {
                calculated(j) = WeightFunction(zvalues(j), InvModel(0));
              }
            double maximum = *std::max_element(calculated.begin(), calculated.end());
            //normalize the predicted data just like the sensitivities
            transform(calculated.begin(), calculated.end(), calculated.begin(),
                [maximum] (double val)
                  { return val/maximum;});
            //calculate the derivatives and the misfit
            for (size_t j = 0; j < zsize; ++j)
              {
                sens(j, 0) = WeightFunction.deriv(zvalues(j), InvModel(0)) / maximum;
                delta(j) = profile(j) - calculated(j);
              }
            // do one step of the inversion
            jif3D::ModelSpaceInversion()(sens, delta, weights, error, 0.0, DeltaModel);
            //and adjust the model (no fancy line search here)
            stepsize = boost::numeric::ublas::norm_2(DeltaModel);
            InvModel -= DeltaModel;
            //Delta model will be used for regularizing the next iteration
            //we don't want this so we set it to zero
            DeltaModel(0) = 0.0;
            ++i;
          }
        for (size_t j = 0; j < zsize; ++j)
          {
            outfile << zvalues(j) << " " << profile(j) << " " << calculated(j)
                << std::endl;
          }
        return InvModel(0);
      }

    /*! Construct a weighting vector that counteracts the decay of the sensitivities of gravity data
     *  to facilitate inversion.
     * @param ZSizes The sizes of the model cells in z-direction in m
     * @param z0 The scaling depth
     * @param WeightVector The resulting weight vector
     * @param WeightFunction The function object that determines the assumes form of the decay
     *
     * The resulting weight vector has only nz elements, where nz is the number of vertical layers
     * the inversion algorithm has to take care to apply these values to the appropriate matrix elements
     */
    void ConstructDepthWeighting(const ThreeDModelBase::t3DModelDim &ZSizes,
        const double z0, rvec &WeightVector, const jif3D::WeightingTerm &WeightFunction)
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
            WeightVector(i) = 1. / WeightFunction(currbottom, z0);

          }
        //normalize
        const double maxel = *std::max_element(WeightVector.begin(), WeightVector.end());
        std::transform(WeightVector.begin(), WeightVector.end(), WeightVector.begin(),
            [maxel] (double val)
              { return val/maxel;});

      }

  }
