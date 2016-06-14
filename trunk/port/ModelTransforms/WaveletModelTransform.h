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
#include "../Global/convert.h"
#include "GeneralModelTransform.h"

namespace jif3D
  {
    //! This class can be used to express the parameters in the wavelet domain and apply regularization there
    /*! This is an experimental transformation class to test whether it works to parametrize the
     * joint inversion model in terms of wavelet parameters. Currently the model has to have
     * a number of cells in each direction that is a power of 2.
     */
    class WaveletModelTransform: public jif3D::GeneralModelTransform
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
      jif3D::rvec Transform(const jif3D::rvec &Vector,
          boost::function1<void, boost::multi_array<double, 3> &> &WaveTrans) const
        {
          if (Vector.size() > Grid.num_elements())
            throw FatalException(
                "In WaveletModelTransform: Model does not fit into grid", __FILE__, __LINE__);

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
          jif3D::rvec result(Vector.size());
          for (size_t j = 0; j < xsize; ++j)
            for (size_t k = 0; k < ysize; ++k)
              {
                for (size_t l = 0; l < zsize; ++l)
                  result(j * (xslicesize) + k * zsize + l) = Grid[j][k][l];
              }
          return result;

        }

    public:
      virtual WaveletModelTransform* clone() const
        {
          return new WaveletModelTransform(*this);
        }
      //! Transform generalized model parameters (wavelet coefficients) to physical parameters
      virtual jif3D::rvec GeneralizedToPhysical(const jif3D::rvec &FullModel) const
        {
          boost::function1<void, boost::multi_array<double, 3> &> func =
              &jif3D::WaveletTransform<boost::multi_array<double, 3> >;
          return Transform(FullModel, func);
        }
      //! Transform physical parameters to generalized model parameters (wavelet coefficients)
      virtual jif3D::rvec PhysicalToGeneralized(const jif3D::rvec &FullModel) const
        {
          boost::function1<void, boost::multi_array<double, 3> &> func =
              &jif3D::InvWaveletTransform<boost::multi_array<double, 3> >;
          return Transform(FullModel, func);
        }
      //! Transform the derivative with respect to physical parameters into the wavelet domain
      virtual jif3D::rvec Derivative(const jif3D::rvec &FullModel,
          const jif3D::rvec &Derivative) const
        {
          return GeneralizedToPhysical(Derivative);
        }
      //! When we construct the model transform, we specify the size of the grid (each size needs to be a power of two)
      WaveletModelTransform(const size_t nx, const size_t ny, const size_t nz) :
          xsize(nx), ysize(ny), zsize(nz)
        {
          if (!IsPowerOfTwo(nx))
            throw jif3D::FatalException(
                "X-dimension of grid is not a power of two: " + jif3D::stringify(nx), __FILE__, __LINE__);
          if (!IsPowerOfTwo(ny))
            throw jif3D::FatalException(
                "Y-dimension of grid is not a power of two: " + jif3D::stringify(ny), __FILE__, __LINE__);
          if (!IsPowerOfTwo(nz))
            throw jif3D::FatalException(
                "Z-dimension of grid is not a power of two: " + jif3D::stringify(nz), __FILE__, __LINE__);
          Grid.resize(boost::extents[nx][ny][nz]);
        }
      };
  }
#endif /* WAVELETMODELTRANSFORM_H_ */
