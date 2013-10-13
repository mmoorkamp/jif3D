//============================================================================
// Name        : WaveletCompressedGravityCalculator.h
// Author      : Dec 12, 2008
// Version     :
// Copyright   : 2008, mmoorkamp
//============================================================================


#ifndef WAVELETCOMPRESSEDGRAVITYCALCULATOR_H_
#define WAVELETCOMPRESSEDGRAVITYCALCULATOR_H_

#include "CachedGravityCalculator.h"

namespace jif3D
  {
    //! Store sensitivity information in a sparse matrix after compression with a 3D wavelet transform and thresholding
    /*! This calculator class performs a 3D wavelet transformation on the row of the sensitivity matrix after
     * the forward calculation for each measurement. It then analyzes the tranformed elements and discards
     * small elements. This decreases the precision of cached forward calculations but greatly reduces the
     * storage requirements.
     *
     * \attention This algorithm can fail spectacularly. There is no threshold value for which a certain
     * precision is guaranteed. In other words even for extremely small threshold values there are some
     * models for which the difference between the exact and approximate result is arbitrarily large (as
     * long as we throw away at least one element). This happens when the transformed and thresholded
     * sensitivities occupy completely different coefficient in the wavelet domain than the model. Use with
     * extreme care !
     */
    class WaveletCompressedGravityCalculator: public jif3D::CachedGravMagCalculator
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
      jif3D::map_mat SparseSens;
      //we try depth weighting to solve some numerical problems
      jif3D::rvec WhiteningVector;
      virtual rvec CalculateNewModel(const ThreeDGravityModel &Model);
      virtual rvec CalculateCachedResult(const ThreeDGravityModel &Model);
      virtual rvec CachedLQDerivative(const ThreeDGravityModel &Model, const rvec &Misfit);
    public:
      //! Get the compressed matrix of sensitivities in the wavelet domain
      const jif3D::map_mat &GetSensitivities() const {return SparseSens;}
      //! Set the desired accuracy for the cached calculations
      /*! More precisely accuracy determines the ratio of the norms of the discarded
       * elements in each row of the sensitivity matrix to the norm of the original
       * row in the sensitivity matrix.
       * @param acc The desired relative accuracy, 0.01 corresponds to 1%
       */
      void SetAccuracy(const double acc){accuracy = acc;}
      //! Set the vector to even out the differences in magnitude between different elements of the sensitivity matrix
      /*! Li and Oldenburg suggest to improve the properties of the compressed matrix
       * to apply the depth weighting before the wavelet transformation. In our
       * experience this does not work.
       * @return A reference to the WhiteningVector
       */
      jif3D::rvec &SetWitheningVector(){return WhiteningVector;}
      virtual void HandleSensitivities(const size_t measindex);
      //! The constructor takes a shared pointer to an implementation object
      WaveletCompressedGravityCalculator(boost::shared_ptr<ThreeDGravMagImplementation> TheImp);
      virtual ~WaveletCompressedGravityCalculator();
      };

  }

#endif /* WAVELETCOMPRESSEDGRAVITYCALCULATOR_H_ */
