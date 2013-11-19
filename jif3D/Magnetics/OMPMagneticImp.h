//============================================================================
// Name        : OMPMagneticImp.h
// Author      : 6 Nov 2013
// Version     : 
// Copyright   : 2013, mm489
//============================================================================

#ifndef OMPMAGNETICIMP_H_
#define OMPMAGNETICIMP_H_

#include <boost/serialization/serialization.hpp>
#include <boost/serialization/base_object.hpp>
#include "../GravMag/ThreeDGravMagImplementation.h"
#include "ThreeDMagneticModel.h"

namespace jif3D
  {

    class OMPMagneticImp: public jif3D::ThreeDGravMagImplementation<ThreeDMagneticModel>
      {
    private:
      double Inclination;
      double Declination;
      double FieldStrength;
      virtual rvec CalcBackground(const size_t measindex, const double xwidth,
          const double ywidth, const double zwidth, const ThreeDMagneticModel &Model,
          rmat &Sensitivities)
      {
          rvec returnvector(ndatapermeas,0.0);
          return returnvector;
      }
      //! Calculate the response of the gridded part
      virtual rvec CalcGridded(const size_t measindex, const ThreeDMagneticModel &Model,
          rmat &Sensitivities);
      //! We have three magnetic field components
      static const size_t ndatapermeas = 3;
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
      OMPMagneticImp(double Inc = 0, double Dec =0, double Fs = 1.0);
      virtual ~OMPMagneticImp();
      };

  } /* namespace jif3D */
#endif /* OMPMAGNETICIMP_H_ */
