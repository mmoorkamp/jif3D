//============================================================================
// Name        : OMPMagneticGradImp.h
// Author      : 17 July 2014
// Version     :
// Copyright   : 2014, zhanjie
//============================================================================

#ifndef OMPMAGNETICGRADIMP_H_
#define OMPMAGNETICGRADIMP_H_

#include "../Global/Serialization.h"
#include "../Global/Jif3DGlobal.h"
#include "../GravMag/ThreeDGravMagImplementation.h"
#include "ThreeDMagneticModel.h"
#include "MagneticData.h"

namespace jif3D
  {

    class J3DEXPORT OMPMagneticGradImp : public jif3D::ThreeDGravMagImplementation<
        MagneticData>
      {
    private:
      double Inclination;
      double Declination;
      double FieldStrength;
      virtual rvec CalcBackground(const size_t measindex, const double xwidth,
          const double ywidth, const double zwidth, const ThreeDMagneticModel &Model,
          const MagneticData &Data,
          rmat &Sensitivities) override
        {
          rvec returnvector(ndatapermeas, 0.0);
          return returnvector;
        }
      //! Calculate the response of the gridded part
      virtual rvec CalcGridded(const size_t measindex, const ThreeDMagneticModel &Model,
          const MagneticData &Data,
          rmat &Sensitivities)  override;
      //! We chose vertical gradient of Z magnetic component from three magnetic field components
      static const size_t ndatapermeas = 1;
      friend class access;
      //! Provide serialization to be able to store objects and, more importantly for simpler MPI parallelization
      template<class Archive>
      void serialize(Archive & ar, const unsigned int version)
        {
          ar & base_object<ThreeDGravMagImplementation<ThreeDMagneticModel>>(*this);
          ar & Inclination;
          ar & Declination;
          ar & FieldStrength;
        }
    public:
      //! How many data do we return before any transformation
      virtual size_t RawDataPerMeasurement() override
        {
          return ndatapermeas;
        }
      OMPMagneticGradImp(double Inc = 0, double Dec = 0, double Fs = 1.0);
      virtual ~OMPMagneticGradImp();
      };

  } /* namespace jif3D */
#endif /* OMPMAGNETICGRADIMP_H_ */
