//============================================================================
// Name        : DepthWeighting.h
// Author      : Sep 18, 2008
// Version     :
// Copyright   : 2008, mmoorkamp
//============================================================================


#ifndef DEPTHWEIGHTING_H_
#define DEPTHWEIGHTING_H_

#include "ThreeDGravityModel.h"

namespace jiba
  {
    /** \addtogroup gravity Gravity forward modelling, display and inversion */
    /* @{ */
    class WeightingTerm
      {
    private:
      double n;
    public:
      double operator()(const double z, const double z0) const;
      double deriv(const double z, const double z0) const;
      double average(const double z1, const double z2, const double z0) const;
      WeightingTerm(const double exponent) :
        n(exponent)
        {
        }
      };

    double FitZ0(const jiba::rvec &SummedSens,
        const jiba::ThreeDGravityModel &Model, const jiba::WeightingTerm &WeightFunction);

    void ConstructDepthWeighting(const ThreeDModelBase::t3DModelDim &XSizes,
        const ThreeDModelBase::t3DModelDim &YSizes,
        const ThreeDModelBase::t3DModelDim &ZSizes, const double z0,
        rvec &WeightVector, const jiba::WeightingTerm &WeightFunction);
    /* @} */
  }
#endif /* DEPTHWEIGHTING_H_ */
