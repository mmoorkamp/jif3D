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
      // the desired accuracy for the cached calculations
      double accuracy;
      //we map the sensitivities on a 3D structure with the same arrangement as the model
      //however each dimension has have a size of power of 2, this structure
      //accomodates this
      boost::multi_array<double, 3> CurrRow;
      //and here are the smallest powers of two that can store the model
      boost::multi_array_types::size_type transformsize[3];
      //the sparse matrix that holds the compressed sensitivities
      jiba::rsparse SparseSens;
      virtual rvec CalculateNewModel(const ThreeDGravityModel &Model);
      virtual rvec CalculateCachedResult(const ThreeDGravityModel &Model);
    public:
      const jiba::rsparse &GetSensitivities() const {return SparseSens;}
      //! Set the desired accuracy for the cached calculations
      /*! More precisely accuracy determines the ratio of the norms of the discarded
       * elements in each row of the sensitivity matrix to the norm of the original
       * row in the sensitivity matrix. Unless the density variations are very large
       * however, this is roughly the accuracy of the result.
       * @param acc The desired relative accuracy, 0.01 corresponds to 1%
       */
      void SetAccuracy(const double acc){accuracy = acc;}
      virtual void HandleSensitivities(const size_t measindex);
      WaveletCompressedGravityCalculator(boost::shared_ptr<ThreeDGravityImplementation> TheImp);
      virtual ~WaveletCompressedGravityCalculator();
      };

  }

#endif /* WAVELETCOMPRESSEDGRAVITYCALCULATOR_H_ */
