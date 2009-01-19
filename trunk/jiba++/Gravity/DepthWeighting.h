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
    /** \addtogroup gravity Gravity forward modeling, display and inversion */
    /* @{ */
    //! This class implements the depth weighting function necessary to counter-act the decay of the sensitivities for gravity data
    /*! Both scalar and tensor gravity data have sensitivity kernels that strongly decay with depth. The result
     * is that in an inversion all structure is concentrated at the top. This class implements functions of the form
     * \f$ (z-z_0)^n \f$ and its derivatives. For scalar data we choose \f$n=-3\f$ and fit \f$ z_0\f$ so we match
     * the decay of the kernel with depth.
     */
    class WeightingTerm
      {
    private:
      double n;
    public:
      //! Given the exponent specified in the constructor calculate the function value for given z and z0
      double operator()(const double z, const double z0) const;
      //! Calculate the derivative of the function for a given z and z0
      double deriv(const double z, const double z0) const;
      //! Calculate the average function value between z1 and z2 for a given z0
      double average(const double z1, const double z2, const double z0) const;
      //! The constructor takes the exponent n
      WeightingTerm(const double exponent) :
        n(exponent)
        {
        }
      };

    //! Given the values of the sensitivity kernel with depth, find z0 that matches the decay
    double FitZ0(const jiba::rvec &SensProfile,
        const ThreeDModelBase::t3DModelDim &ZSizes, const jiba::WeightingTerm &WeightFunction);
    //! Given a z0 and the model geometry, construct a vector of weights
    void ConstructDepthWeighting(
        const ThreeDModelBase::t3DModelDim &ZSizes, const double z0,
        rvec &WeightVector, const jiba::WeightingTerm &WeightFunction);
    //! Extract sensitivities for a site that is closest to the middle of the modeling domain
    void ExtractMiddleSens(const jiba::ThreeDGravityModel &Model,
        const jiba::rmat &Sensitivities, const size_t MeasPerPos, jiba::rvec &SensProfile);
    /* @} */
  }
#endif /* DEPTHWEIGHTING_H_ */
