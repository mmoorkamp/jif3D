//============================================================================
// Name        : WaveletModelTransform.h
// Author      : Apr 4, 2011
// Version     : 
// Copyright   : 2011, mmoorkamp
//============================================================================


#ifndef WAVELETMODELTRANSFORM_H_
#define WAVELETMODELTRANSFORM_H_

#include <boost/multi_array.hpp>
#include <boost/function.hpp>
#include "../Global/FatalException.h"
#include "../Global/Wavelet.h"
#include "../Global/NumUtil.h"
#include "ModelDistributor.h"

namespace jiba
  {
    //! This class can be used to express the parameters in the wavelet domain and apply regularization there
    /*! This is an experimental transformation class to test whether it works to parametrize the
     * joint inversion model in terms of wavelet parameters. Currently the model has to have
     * a number of cells in each direction that is a power of 2.
     */
    class WaveletModelTransform: public jiba::GeneralModelTransform
      {
    private:
      //! The size of the grid in x-direction
      const size_t xsize;
      //! The size of the grid in y-direction
      const size_t ysize;
      //! The size of the grid in z-direction
      const size_t zsize;
      //! We need to transform the vector of model parameters to a 3D grid with correct geometry, this member is used to create this grid
      mutable boost::multi_array<double, 3> Grid;
      //! Copy the model vector to the 3D grid and perform either a forward or inverse wavelet transform, then copy back to a vector
      jiba::rvec Transform(const jiba::rvec &Vector, boost::function1<void,
          boost::multi_array<double, 3> > &WaveTrans) const
        {
          if (Vector.size() > Grid.num_elements())
            throw FatalException(
                "In WaveletModelTransform: Model does not fit into grid");

          //FullModel and Grid have different sizes because we are padding
          // Grid to the nearest power of two, but we have to
          //consider the spatial structure of the underlying model, so we have
          //to use loops for copying and NOT std::copy
          const size_t xslicesize = ysize * zsize;
          std::fill_n(Grid.origin(), Grid.num_elements(), 0.0);
          for (size_t j = 0; j < xsize; ++j)
            for (size_t k = 0; k < ysize; ++k)
              {
                for (size_t l = 0; l < zsize; ++l)
                  Grid[j][k][l] = Vector(j * (xslicesize) + k * zsize + l);
              }
          WaveTrans(Grid);
          jiba::rvec result(Vector.size());
          for (size_t j = 0; j < xsize; ++j)
            for (size_t k = 0; k < ysize; ++k)
              {
                for (size_t l = 0; l < zsize; ++l)
                  result(j * (xslicesize) + k * zsize + l) = Grid[j][k][l];
              }
          return result;

        }

    public:
      //! Transform generalized model parameters (wavelet coefficients) to physical parameters
      virtual jiba::rvec GeneralizedToPhysical(const jiba::rvec &FullModel) const
        {
          boost::function1<void, boost::multi_array<double, 3> > func =
              &jiba::InvWaveletTransform<boost::multi_array<double, 3> >;
          return Transform(FullModel, func);
        }
      //! Transform physical parameters to generalized model parameters (wavelet coefficients)
      virtual jiba::rvec PhysicalToGeneralized(const jiba::rvec &FullModel) const
        {
          boost::function1<void, boost::multi_array<double, 3> > func =
              &jiba::WaveletTransform<boost::multi_array<double, 3> >;
          return Transform(FullModel, func);
        }
      //! Transform the derivative with respect to physical parameters into the wavelet domain
      virtual jiba::rvec Derivative(const jiba::rvec &FullModel,
          const jiba::rvec &Derivative) const
        {
          return GeneralizedToPhysical(Derivative);
        }
      //! When we construct the model transform, we specify the size of the grid (each size needs to be a power of two)
      WaveletModelTransform(const size_t nx, const size_t ny, const size_t nz) :
        xsize(nx), ysize(ny), zsize(nz)
        {
          if (!IsPowerOfTwo(nx))
            throw jiba::FatalException(
                "X-dimension of grid is not a power of two.");
          if (!IsPowerOfTwo(ny))
            throw jiba::FatalException(
                "Y-dimension of grid is not a power of two.");
          if (!IsPowerOfTwo(nz))
            throw jiba::FatalException(
                "Z-dimension of grid is not a power of two.");
          Grid.resize(boost::extents[nx][ny][nz]);
        }
      };
  }
#endif /* WAVELETMODELTRANSFORM_H_ */
