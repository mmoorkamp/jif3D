//============================================================================
// Name        : DepthWeigthing.cpp
// Author      : Sep 18, 2008
// Version     :
// Copyright   : 2008, mmoorkamp
//============================================================================

#include "DepthWeighting.h"
#include "../Inversion/LinearInversion.h"
#include <fstream>
#include <iostream>

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

    double FitZ0(const jiba::rvec &SummedSens,
        const jiba::ThreeDGravityModel &Model, const jiba::WeightingTerm &WeightFunction)
      {
        const size_t xsize = Model.GetXCellSizes().shape()[0];
        const size_t ysize = Model.GetYCellSizes().shape()[0];
        const size_t zsize = Model.GetZCellSizes().shape()[0];
        //the index of the center horizontal cell
        size_t startindex = zsize * (ysize -1) * xsize /2;
        jiba::rvec sensprofile(zsize), zvalues(zsize);
        partial_sum(Model.GetZCellSizes().begin(), Model.GetZCellSizes().end(),
            zvalues.begin());
        copy(SummedSens.begin() + startindex, SummedSens.begin() + startindex
            + zsize, sensprofile.begin());
        transform(sensprofile.begin(), sensprofile.end(), sensprofile.begin(),
            boost::bind(std::divides<double>(), _1, *std::max_element(
                    sensprofile.begin(), sensprofile.end())));
        transform(sensprofile.begin(), sensprofile.end(), sensprofile.begin(),std::abs<double>);
        double startz = zvalues(zsize - 1) * 1.1;

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
        while (stepsize> 0.1 && i < iterations)
          {
            for (size_t j = 0; j < zsize; ++j)
              {
                calculated(j) = WeightFunction(zvalues(j),InvModel(0));
                delta(j) = sensprofile(j) - calculated(j);
              }
            double maximum = *std::max_element(calculated.begin(),
                calculated.end());
            transform(calculated.begin(), calculated.end(), calculated.begin(),
                boost::bind(std::divides<double>(), _1, maximum));
            for (size_t j = 0; j < zsize; ++j)
              {
                sens(j, 0) = WeightFunction.deriv(zvalues(j), InvModel(0)) / maximum;
                delta(j) = sensprofile(j) - calculated(j);
              }
            for (size_t j = 0; j < zsize; ++j)
              {
                outfile << zvalues(j) << " " << sensprofile(j) << " "
                << calculated(j) << std::endl;
              }

            jiba::DataSpaceInversion()(sens, delta, weights, error, evalthresh,
                0.0, DeltaModel);

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
