/*
 * CudaMagneticSusceptibilityImp.h
 *
 *  Created on: Jul 19, 2022
 *      Author: max
 */

#ifndef MAGNETICS_CUDAMAGNETICSUSCEPTIBILITYIMP_H_
#define MAGNETICS_CUDAMAGNETICSUSCEPTIBILITYIMP_H_

#include "OMPMagneticSusceptibilityImp.h"

namespace jif3D
  {

    class CudaMagneticSusceptibilityImp: public OMPMagneticSusceptibilityImp
      {
    private:
      /*! These pointers hold the memory on the graphics card as allocated
       * before the calculation. We make them class variables so that we
       * can use them in different parts of the class and only have to
       * do the allocation once for all measurements
       * */
      double *d_xcoord, *d_ycoord, *d_zcoord;
      //! Pointers for the cell size structures on the GPU
      double *d_xsize, *d_ysize, *d_zsize;
      //! The pointer to the array of sensitivies on the GPU for the current tensor element
      double *d_result;
      //! We need a raw double pointer to store the sensitivities for  the current measurements, see CalcGridded
      double *currsens;
      //! The number of elements in one row of the sensitivity matrix
      size_t currsenssize;
      //! The size of a CUDA execution block
      size_t blocksize;
    public:
      virtual rvec Calculate(const ThreeDSusceptibilityModel &Model, const TotalFieldMagneticData &Data,
          ThreeDGravMagCalculator<TotalFieldMagneticData> &Calculator) override;
      virtual rvec CalcGridded(const size_t measindex,
          const ThreeDSusceptibilityModel &Model, const TotalFieldMagneticData &Data,
          rmat &Sensitivities) override;
      CudaMagneticSusceptibilityImp(double Inc = 0, double Dec = 0, double Fs = 1.0);
      virtual ~CudaMagneticSusceptibilityImp();
      };

  } /* namespace jif3D */

#endif /* MAGNETICS_CUDAMAGNETICSUSCEPTIBILITYIMP_H_ */
