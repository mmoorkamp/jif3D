//============================================================================
// Name        : WaveletCompressedGravityCalculator.h
// Author      : Dec 12, 2008
// Version     :
// Copyright   : 2008, mmoorkamp
//============================================================================


#ifndef WAVELETCOMPRESSEDGRAVITYCALCULATOR_H_
#define WAVELETCOMPRESSEDGRAVITYCALCULATOR_H_

#include "CachedGravityCalculator.h"

namespace jiba
  {

    class WaveletCompressedGravityCalculator: public jiba::CachedGravityCalculator
      {
    private:
      size_t nmeas ;
      size_t ngrid ;
      size_t nbglayers;
      size_t nmod;
      size_t xsize ;
      size_t ysize ;
      size_t zsize ;
      boost::multi_array<double, 3> CurrRow;
      boost::multi_array_types::size_type transformsize[3];
      boost::numeric::ublas::mapped_matrix<double,
          boost::numeric::ublas::column_major> SparseSens;
      virtual rvec CalculateNewModel(const ThreeDGravityModel &Model);
      virtual rvec CalculateCachedResult(const ThreeDGravityModel &Model);
    public:
      virtual void HandleSensitivities(const size_t measindex);

      WaveletCompressedGravityCalculator(ThreeDGravityImplementation &TheImp);
      virtual ~WaveletCompressedGravityCalculator();
      };

  }

#endif /* WAVELETCOMPRESSEDGRAVITYCALCULATOR_H_ */
