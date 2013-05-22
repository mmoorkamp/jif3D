//============================================================================
// Name        : TensorCudaGravityImp.h
// Author      : Feb 11, 2009
// Version     :
// Copyright   : 2009, mmoorkamp
//============================================================================

#ifndef TENSORCUDAGRAVITYIMP_H_
#define TENSORCUDAGRAVITYIMP_H_

#include <boost/serialization/serialization.hpp>
#include <boost/serialization/base_object.hpp>
#include "ThreeDGravityImplementation.h"

namespace jif3D
  {

    /** \addtogroup gravity Gravity forward modeling, display and inversion */
    /* @{ */
    //! Calculate a FTG gravity response using Cuda
    /*! This class is the tensorial counterpart to ScalarCudaGravityImp.
     * It calculates the 9 elements of second derivatives of the gravitational potential.
     * It only implements the calculation of the background and the gridded part.
     * The assembly of the two parts is performed by the base class ThreeDGravityImplementation.
     *
     * This implementation class uses Nvidia's CUDA API to perform
     * the forward calculation on a Nvidia graphics card. It needs
     * a card with compute capability 1.3 or more to perform double
     * length floating point calculations.
     *
     * This class cannot be copied as this would make a mess with management
     * of cuda resources.
     */
    class TensorCudaGravityImp: public jif3D::ThreeDGravityImplementation
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
      friend class boost::serialization::access;
      //! Provide serialization to be able to store objects and, more importantly for simpler MPI parallelization
      template<class Archive>
      void serialize(Archive & ar, const unsigned int version)
        {
          ar & boost::serialization::base_object<ThreeDGravityImplementation>(*this);
          ar & d_xcoord;
          ar & d_ycoord;
          ar & d_zcoord;
          ar & d_xsize;
          ar & d_ysize;
          ar & d_zsize;
          ar & d_result;
          ar & currsens;
          ar & currsenssize;
          ar & blocksize;
        }
      //! This is a tensor calculation for any outside call we return 9 data per measurement
      static const size_t ndatapermeas = 9;
      //! Calculate the response of the background, currently this is done on the CPU
      virtual rvec CalcBackground(const size_t measindex, const double xwidth,
          const double ywidth, const double zwidth, const ThreeDGravityModel &Model,
          rmat &Sensitivities);
      //! Calculate the response of the gridded part, this is done on the GPU with CUDA
      virtual rvec CalcGridded(const size_t measindex, const ThreeDGravityModel &Model,
          rmat &Sensitivities);
      //! This class cannot be copied, so copy constructor and assignment are private
      jif3D::TensorCudaGravityImp &operator=(const jif3D::TensorCudaGravityImp&);
      //! This class cannot be copied, so copy constructor and assignment are private
      TensorCudaGravityImp(const jif3D::TensorCudaGravityImp&);
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
          ThreeDGravityCalculator &Calculator);
      TensorCudaGravityImp();
      virtual ~TensorCudaGravityImp();
      };

  }

#endif /* TENSORCUDAGRAVITYIMP_H_ */
