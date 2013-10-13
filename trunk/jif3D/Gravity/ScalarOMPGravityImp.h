//============================================================================
// Name        : ScalarOMPGravityImp.h
// Author      : Dec 10, 2008
// Version     :
// Copyright   : 2008, mmoorkamp
//============================================================================


#ifndef SCALAROMPGRAVITYIMP_H_
#define SCALAROMPGRAVITYIMP_H_

#include <boost/serialization/serialization.hpp>
#include <boost/serialization/base_object.hpp>
#include "../GravMag/ThreeDGravMagImplementation.h"

namespace jif3D
  {
    /** \addtogroup gravity Gravity forward modeling, display and inversion */
    /* @{ */
    //! Calculate a scalar gravity response using OpenMP parallelization
    /*! Calculation of scalar gravity acceleration using parallelization with openmp.
     * This is the standard implementation and as openmp uses pragma directives can
     * also be used in serial with a compiler that does not understand openmp.
     */
    class ScalarOMPGravityImp: public jif3D::ThreeDGravMagImplementation
      {
    private:
      //! We calculate scalar data, so we have one data point per site
      static const size_t ndatapermeas = 1;
      //! Implement the calculation of the background response
      virtual rvec CalcBackground(const size_t measindex, const double xwidth, const double ywidth,
          const double zwidth, const ThreeDGravityModel &Model,
          rmat &Sensitivities);
      //! Calculate the response of the gridded part
      virtual rvec CalcGridded(const size_t measindex, const ThreeDGravityModel &Model,
          rmat &Sensitivities);
      friend class boost::serialization::access;
      //! Provide serialization to be able to store objects and, more importantly for simpler MPI parallelization
      template<class Archive>
      void serialize(Archive & ar, const unsigned int version)
        {
          ar & boost::serialization::base_object<ThreeDGravMagImplementation>(*this);
        }
    public:
      //! How many data do we return before any transformation
      virtual size_t RawDataPerMeasurement()
        {
          return ndatapermeas;
        }
      ScalarOMPGravityImp();
      virtual ~ScalarOMPGravityImp();
      };
  /* @} */
  }

#endif /* SCALAROMPGRAVITYIMP_H_ */
