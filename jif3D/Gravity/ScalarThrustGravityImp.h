//============================================================================
// Name        : ScalarThrustGravityImp.h
// Author      : 15 Dec 2016
// Version     : 
// Copyright   : 2016, mm489
//============================================================================

#ifndef GRAVITY_SCALARTHRUSTGRAVITYIMP_H_
#define GRAVITY_SCALARTHRUSTGRAVITYIMP_H_

#include "../Global/Serialization.h"
#include "../Global/Jif3DGlobal.h"
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>

namespace jif3D
  {

    class ScalarThrustGravityImp
      {
    private:
      //! This is a scalar calculation so we get one value per measurement
      static const size_t ndatapermeas = 1;
      thrust::device_vector<double> d_xcoord, d_ycoord, d_zcoord;
      thrust::device_vector<double> d_xsize, d_ysize, d_zsize;
      thrust::device_vector<double> d_result;
      //! Calculate the response of the background, currently this is done on the CPU
      virtual rvec CalcBackground(const size_t measindex, const double xwidth,
          const double ywidth, const double zwidth,
          const ThreeDGravityModel &Model, rmat &Sensitivities);
      //! Calculate the response of the gridded part, this is done on the GPU with CUDA
      virtual rvec CalcGridded(const size_t measindex,
          const ThreeDGravityModel &Model, rmat &Sensitivities);
    public:
      //! How many data do we return before any transformation
      virtual size_t RawDataPerMeasurement()
        {
          return ndatapermeas;
        }
      //! Set the size for a CUDA execution block, the default value is 128
      void SetCUDABlockSize(const size_t s)
        {
          blocksize = s;
        }
      //! We reimplement the Calculate method to accommodate some specific CUDA issues
      virtual rvec Calculate(const ThreeDGravityModel &Model,
          ThreeDGravMagCalculator<ThreeDGravityModel> &Calculator);
      ScalarThrustGravityImp();
      virtual ~ScalarThrustGravityImp();
      };

  } /* namespace jif3D */

#endif /* GRAVITY_SCALARTHRUSTGRAVITYIMP_H_ */
