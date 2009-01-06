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
    //! Store sensitivity information in a sparse matrix after compression with a 3D wavelet transform and thresholding
    /*! This calculator class performs a 3D wavelet transformation on the row of the sensitivity matrix after
     * the forward calculation for each measurement. It then analyzes the tranformed elements and discards
     * small elements. This decreases the precision of cached forward calculations but greatly reduces the
     * storage requirements.
     */
    class WaveletCompressedGravityCalculator: public jiba::CachedGravityCalculator
      {
    private:
      //we store some array sizes when we calculate a new model, so we do not
      //have to recalculate them when we do a cached calculation
      size_t nmeas ;
      size_t ngrid ;
      size_t nbglayers;
      size_t nmod;
      size_t xsize ;
      size_t ysize ;
      size_t zsize ;
      //we map the sensitivities on a 3D structure with the same arrangement as the model
      //however each dimension has have a size of power of 2, this structure
      //accomodates this
      boost::multi_array<double, 3> CurrRow;
      //and here are the smallest powers of two that can store the model
      boost::multi_array_types::size_type transformsize[3];
      //the sparse matrix that holds the compressed sensitivities
      boost::numeric::ublas::mapped_matrix<double,
          boost::numeric::ublas::column_major> SparseSens;
      virtual rvec CalculateNewModel(const ThreeDGravityModel &Model);
      virtual rvec CalculateCachedResult(const ThreeDGravityModel &Model);
    public:
      virtual void HandleSensitivities(const size_t measindex);
      WaveletCompressedGravityCalculator(boost::shared_ptr<ThreeDGravityImplementation> TheImp);
      virtual ~WaveletCompressedGravityCalculator();
      };

  }

#endif /* WAVELETCOMPRESSEDGRAVITYCALCULATOR_H_ */
