/*
 * ScalarCudaGravityImp.h
 *
 *  Created on: Dec 10, 2008
 *      Author: mmoorkamp
 */

#ifndef SCALARCUDAGRAVITYIMP_H_
#define SCALARCUDAGRAVITYIMP_H_

#include "ThreeDGravityImplementation.h"

namespace jiba
  {
    /** \addtogroup gravity Gravity forward modelling, display and inversion */
    /* @{ */
    //! Calculate a scalar gravity response using Nvidia's CUDA API
    /*! This implementation class uses Nvidia's CUDA api to perform
     * the forward calculation on a Nvidia graphics card. It needs
     * a card with compute capability 1.3 or more to perform double
     * length floating point calculations
     */
    class ScalarCudaGravityImp: public jiba::ThreeDGravityImplementation
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
      // This is a scalar calculation so we get one value per measurement
      static const size_t ndatapermeas = 1;
      //! Calculate the response of the background, currently this is done on the CPU
      virtual rvec CalcBackground(const size_t measindex, const double xwidth, const double ywidth,
          const double zwidth, const ThreeDGravityModel &Model,
          rmat &Sensitivities);
      //! Calculate the response of the gridded part, this is done on the GPU with CUDA
      virtual rvec CalcGridded(const size_t measindex, const ThreeDGravityModel &Model,
          rmat &Sensitivities);
    public:
      virtual size_t GetDataPerMeasurement()
        {
          return ndatapermeas;
        }
      //! We reimplement the Calculate method to accommodate some specific CUDA issues
      virtual rvec Calculate(const ThreeDGravityModel &Model,
          ThreeDGravityCalculator &Calculator);
      ScalarCudaGravityImp();
      virtual ~ScalarCudaGravityImp();
      };
  /* @} */
  }

#endif /* SCALARCUDAGRAVITYIMP_H_ */
