/*
 * OMPMagnetizationImp.h
 *
 *  Created on: May 19, 2022
 *      Author: max
 */

#ifndef MAGNETICS_OMPMAGNETIZATIONIMP_H_
#define MAGNETICS_OMPMAGNETIZATIONIMP_H_

#include "../Global/Serialization.h"
#include "../Global/Jif3DGlobal.h"
#include "../GravMag/ThreeDGravMagImplementation.h"
#include "ThreeDMagnetizationModel.h"
#include "ThreeComponentMagneticData.h"

namespace jif3D
  {

    class OMPMagnetizationImp: public ThreeDGravMagImplementation<
        ThreeComponentMagneticData>
      {
    private:
      virtual rvec CalcBackground(const size_t measindex, const double xwidth,
          const double ywidth, const double zwidth,
          const ThreeDMagnetizationModel &Model, const ThreeComponentMagneticData &Data,
          rmat &Sensitivities) override;
      //! Calculate the response of the gridded part
      virtual rvec CalcGridded(const size_t measindex,
          const ThreeDMagnetizationModel &Model, const ThreeComponentMagneticData &Data,
          rmat &Sensitivities) override;
      //! We have three magnetic field components
      static const size_t ndatapermeas = 3;
      friend class access;
      //! Provide serialization to be able to store objects and, more importantly for simpler MPI parallelization
      template<class Archive>
      void serialize(Archive &ar, const unsigned int version)
        {
          ar & base_object<ThreeDGravMagImplementation>(*this);

        }
    public: //! How many data do we return before any transformation
      virtual size_t RawDataPerMeasurement() override
        {
          return ndatapermeas;
        }
      OMPMagnetizationImp();
      virtual ~OMPMagnetizationImp();
      };

  } /* namespace jif3D */

#endif /* MAGNETICS_OMPMAGNETIZATIONIMP_H_ */
