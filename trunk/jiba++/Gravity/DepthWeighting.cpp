//============================================================================
// Name        : DepthWeigthing.cpp
// Author      : Sep 18, 2008
// Version     :
// Copyright   : 2008, mmoorkamp
//============================================================================

#include <fstream>
#include <iostream>
#include <boost/bind.hpp>

#include "DepthWeighting.h"
#include "../Inversion/LinearInversion.h"

namespace jiba
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
     */
    double WeightingTerm::deriv(const double z, const double z0) const
      {
        return n * pow(z + z0, n - 1);
      }
    /*! Calculate the average value of f in the intervall [z1,z2]
     * @param z1 The lower bound of the intervall
     * @param z2 The upper bound of the intervall
     * @param z0 The parametric depth value
     * @return The average value of f
     */
    double WeightingTerm::average(const double z1, const double z2,
        const double z0) const
      {
        return 1. / ((z2 - z1) * (n + 1)) * (pow(z2 + z0, n + 1) - pow(z1 + z0,
            n + 1));
      }
    /*!
     * @param SensProfile A vector with nz elements that contains the sensitivities (for a single measurement) for each cell with depth
     * @param ZSizes The size of each cell in z-direction
     * @param WeightFunction The function object that gives the form of the weighting terms and it derivatives, e.g. an instance of the class WeightingTerm
     * @return An estimate of z0 for which the weighting function matches the sensitivity decay most closely
     */
    double FitZ0(const jiba::rvec &SensProfile,
        const ThreeDModelBase::t3DModelDim &ZSizes,
        const jiba::WeightingTerm &WeightFunction)
      {

        const size_t zsize = SensProfile.size();
        jiba::rvec profile(SensProfile), zvalues(zsize);
        //we need the depth to the interfaces
        partial_sum(ZSizes.begin(), ZSizes.end(), zvalues.begin());
        //the sensitivities can be negative, so we examine the absolute value
        transform(profile.begin(), profile.end(), profile.begin(), std::abs<
            double>);

        //then we normalize so the largest element (usually the top cell) has value 1
        transform(profile.begin(), profile.end(), profile.begin(), boost::bind(
            std::divides<double>(), _1, *std::max_element(profile.begin(),
                profile.end())));
        // we guess a value for z0, we take it slightly larger than the bottom of the model
        //it shouldn't be on a cell boundary to avoid numerical problems
        double startz = zvalues(zsize - 1) * 1.1;

        //setup the "inversion"
        const size_t ndata = zvalues.size();
        jiba::rmat sens(zsize, 1);
        jiba::rvec weights(1);
        jiba::rvec InvModel(1), DeltaModel(1);
        weights(0) = 1.0;
        InvModel(0) = startz;
        DeltaModel(0) = 0.0;
        const size_t iterations = 100;
        jiba::rvec error(ndata), delta(ndata), calculated(ndata);
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
                calculated( j) = WeightFunction(zvalues(j), InvModel(0));
              }
            double maximum = *std::max_element(calculated.begin(),
                calculated.end());
            //normalize the predicted data just like the sensitivities
            transform(calculated.begin(), calculated.end(), calculated.begin(),
                boost::bind(std::divides<double>(), _1, maximum));
            //calculate the derivatives and the misfit
            for (size_t j = 0; j < zsize; ++j)
              {
                sens(j, 0) = WeightFunction.deriv(zvalues(j), InvModel(0))
                    / maximum;
                delta( j) = profile(j) - calculated(j);
              }
            // do one step of the inversion
            jiba::ModelSpaceInversion()(sens, delta, weights, error, 0.0, DeltaModel);
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
        const double z0, rvec &WeightVector,
        const jiba::WeightingTerm &WeightFunction)
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
            WeightVector( i) = 1. / WeightFunction(currbottom, z0);

          }
        //normalize
        std::transform(WeightVector.begin(), WeightVector.end(),
            WeightVector.begin(), boost::bind(std::divides<double>(), _1,
                *std::max_element(WeightVector.begin(), WeightVector.end())));

      }
    /*! For depth weighting we need the sensitivities below a site to match with our weighting function.
     * We somewhat arbitrarily chose the site closest to the middle, as here the effects of the finite modeling domain
     * should be smallest. The function extracts the row from the matrix that corresponds to this location
     * @param Model The model object, needed for the grid information
     * @param Sensitivities The sensitivity matrix that we want to extract the profile from
     * @param MeasPerPos How many measurements per site,e.g. 1 for only scalar and 9 for only FTG
     * @param SensProfile The depth profile of the sensitivity below the site
     */
    void ExtractMiddleSens(const jiba::ThreeDGravityModel &Model,
        const jiba::rmat &Sensitivities, const size_t MeasPerPos,
        jiba::rvec &SensProfile)
      {
    	const size_t nmeas = Model.GetMeasPosX().size();
        const double midx =
            Model.GetXCoordinates()[Model.GetXCoordinates().size() - 1] / 2.0;
        const double midy =
            Model.GetYCoordinates()[Model.GetYCoordinates().size() - 1] / 2.0;


        jiba::rvec distances(nmeas);
        for (size_t i = 0; i < nmeas; ++i)
          {
            distances( i) = sqrt(pow(Model.GetMeasPosX()[i] - midx, 2) + pow(
                Model.GetMeasPosY()[i] - midy, 2));
          }
        const size_t midindex = distance(distances.begin(), std::min_element(
            distances.begin(), distances.end()));
        boost::array<jiba::ThreeDModelBase::t3DModelData::index, 3> modelindex(
            Model.FindAssociatedIndices(Model.GetMeasPosX()[midindex],
                Model.GetMeasPosY()[midindex], 0.0));
        //we store the sensitivities for the background at the end of the matrix
        //so we can ignore it here
        boost::numeric::ublas::matrix_row<const jiba::rmat> MiddleSens(Sensitivities, midindex * MeasPerPos);

        const size_t ysize = Model.GetDensities().shape()[1];
        const size_t zsize = Model.GetDensities().shape()[2];

        SensProfile.resize(zsize);
        //the same here, if we operate on the first ngrid elements
        //the background does not matter
        const size_t startindex = (zsize * ysize) * modelindex[0] + zsize
            * modelindex[1];
        std::copy(MiddleSens.begin() + startindex, MiddleSens.begin()
            + startindex + zsize, SensProfile.begin());

      }
  }
