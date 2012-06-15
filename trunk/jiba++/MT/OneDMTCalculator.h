//============================================================================
// Name        : OneDMTCalculator.h
// Author      : 13 Jun 2012
// Version     : 
// Copyright   : 2012, mm489
//============================================================================

#ifndef ONEDMTCALCULATOR_H_
#define ONEDMTCALCULATOR_H_

#include "../Global/VecMat.h"
#include "X3DModel.h"

namespace jiba
  {
    /** \addtogroup mtmodelling Forward modelling of magnetotelluric data */
    /* @{ */
    class OneDMTCalculator
      {
    public:
      //! This type definition is necessary so that ThreeDModelObjective can correctly deduce the native type for a model object for this class
      typedef X3DModel ModelType;
    private:
      jiba::cmat Gamma;
      jiba::cmat gammakj;
      jiba::cmat gammaj;
      jiba::cvec Z;
    public:
      const jiba::cmat &GetGamma()
        {
          return Gamma;
        }
      rvec Calculate(const ModelType &Model);
      rvec LQDerivative(const ModelType &Model, const rvec &Misfit);
      OneDMTCalculator();
      virtual ~OneDMTCalculator();
      };
  /* @} */
  } /* namespace jiba */
#endif /* ONEDMTCALCULATOR_H_ */
