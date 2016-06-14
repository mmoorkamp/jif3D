//============================================================================
// Name        : ScalarCudaGravityImp.h
// Author      : Dec 10, 2008
// Version     :
// Copyright   : 2009, mmoorkamp
//============================================================================

#ifndef SCALARCUDAGRAVITYIMP_H_
#define SCALARCUDAGRAVITYIMP_H_

#include "../Global/Serialization.h"
#include "../Global/Jif3DGlobal.h"
#include "../Gravity/ThreeDGravityModel.h"
#include "../GravMag/ThreeDGravMagImplementation.h"

namespace jif3D
  {
    /** \addtogroup gravity Gravity forward modeling, display and inversion */
    /* @{ */
    //! Calculate a scalar gravity response using Nvidia's CUDA API
    /*! This implementation class uses Nvidia's CUDA API to perform
     * the forward calculation on a Nvidia graphics card. It needs
     * a card with compute capability 1.3 or more to perform double
     * length floating point calculations.
     *
     * This class cannot be copied as this would make a mess with management
     * of cuda resources.
     */
    class J3DEXPORT ScalarCudaGravityImp: public jif3D::ThreeDGravMagImplementation<ThreeDGravityModel>
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
      //! The size of the current sensitivity row
      size_t currsenssize;
      //! The size of a CUDA execution block
      size_t blocksize;
      //! This is a scalar calculation so we get one value per measurement
      static const size_t ndatapermeas = 1;
      friend class access;
      //! Provide serialization to be able to store objects and, more importantly for simpler MPI parallelization
      template<class Archive>
      void serialize(Archive & ar, const unsigned int version)
        {
          ar & boost::serialization::base_object<ThreeDGravMagImplementation>(*this);
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
      //! Calculate the response of the background, currently this is done on the CPU
      virtual rvec CalcBackground(const size_t measindex, const double xwidth,
          const double ywidth, const double zwidth,
          const ThreeDGravityModel &Model, rmat &Sensitivities);
      //! Calculate the response of the gridded part, this is done on the GPU with CUDA
      virtual rvec CalcGridded(const size_t measindex,
          const ThreeDGravityModel &Model, rmat &Sensitivities);
      //! This class cannot be copied, so copy constructor and assignment are private
      jif3D::ScalarCudaGravityImp &operator=(const jif3D::ScalarCudaGravityImp&);
      ScalarCudaGravityImp(const jif3D::ScalarCudaGravityImp&);
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
      ScalarCudaGravityImp();
      virtual ~ScalarCudaGravityImp();
      };
  /* @} */
  }

#endif /* SCALARCUDAGRAVITYIMP_H_ */
