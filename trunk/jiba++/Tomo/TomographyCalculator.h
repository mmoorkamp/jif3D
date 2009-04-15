//============================================================================
// Name        : TomographyCalculator.h
// Author      : Apr 14, 2009
// Version     :
// Copyright   : 2009, mmoorkamp
//============================================================================


#ifndef TOMOGRAPHYCALCULATOR_H_
#define TOMOGRAPHYCALCULATOR_H_

#include "ThreeDSeismicModel.h"
#include "../Global/VecMat.h"
#include "modeling_seismic.h"

namespace jiba
  {

    class TomographyCalculator
      {
      private:
        jiba::GEOMETRY geo;
        jiba::GRID_STRUCT grid;
        jiba::DATA_STRUCT data;
        jiba::RP_STRUCT *raypath;
        void Allocate(const size_t ngrid,const size_t ndata,const size_t npos);
    public:
      TomographyCalculator();
      virtual ~TomographyCalculator();
      rvec Calculate(const ThreeDSeismicModel &Model);
      rvec LQDerivative(const ThreeDSeismicModel &Model, const rvec &Misfit);
      };

  }

#endif /* TOMOGRAPHYCALCULATOR_H_ */
