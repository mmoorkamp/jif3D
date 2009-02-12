//============================================================================
// Name        : TensorCudaGravityImp.h
// Author      : Feb 11, 2009
// Version     :
// Copyright   : 2009, mmoorkamp
//============================================================================


#ifndef TENSORCUDAGRAVITYIMP_H_
#define TENSORCUDAGRAVITYIMP_H_

#include "ThreeDGravityImplementation.h"

namespace jiba
  {

    class TensorCudaGravityImp: public jiba::ThreeDGravityImplementation
      {
    private:
      // These pointers hold the memory on the graphics card as allocated
      // before the calculation. We make them class variables so that we
      // can use them in different parts of the program and only have to
      // do the allocation once for all measurements
      double *d_xcoord, *d_ycoord, *d_zcoord;
      double *d_xsize, *d_ysize, *d_zsize;
      double *d_result;
      // we need a raw double pointer to store the sensitivities for
      // the current measurements, see CalcGridded
      double *currsens;
      size_t currsenssize;
      // This is a tensor calculation for any outside call we return 9 data per measurement
      static const size_t ndatapermeas = 9;
      //! Calculate the response of the background, currently this is done on the CPU
      virtual rvec CalcBackground(const size_t measindex, const double xwidth,
          const double ywidth, const double zwidth,
          const ThreeDGravityModel &Model, rmat &Sensitivities);
      //! Calculate the response of the gridded part, this is done on the GPU with CUDA
      virtual rvec CalcGridded(const size_t measindex,
          const ThreeDGravityModel &Model, rmat &Sensitivities);
    public:
      virtual size_t GetDataPerMeasurement()
        {
          return ndatapermeas;
        }
      //! We reimplement the Calculate method to accommodate some specific CUDA issues
      virtual rvec Calculate(const ThreeDGravityModel &Model,
          ThreeDGravityCalculator &Calculator);
      TensorCudaGravityImp();
      virtual ~TensorCudaGravityImp();
      };

  }

#endif /* TENSORCUDAGRAVITYIMP_H_ */
