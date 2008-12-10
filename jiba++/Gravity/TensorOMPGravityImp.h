//============================================================================
// Name        : TensorOMPGravityImp.h
// Author      : Dec 10, 2008
// Version     :
// Copyright   : 2008, mmoorkamp
//============================================================================


#ifndef TENSOROMPGRAVITYIMP_H_
#define TENSOROMPGRAVITYIMP_H_

#include "ThreeDGravityImplementation.h"

namespace jiba
  {

    class TensorOMPGravityImp: public jiba::ThreeDGravityImplementation
      {
    private:
      static const size_t ndatapermeas = 9;
    public:
      virtual size_t GetDataPerMeasurement()
        {
          return ndatapermeas;
        }
      virtual rvec CalcBackground(const double xmeas, const double ymeas,
          const double zmeas, const double xwidth, const double ywidth,
          const double zwidth, const ThreeDGravityModel &Model,
          rmat &Sensitivities);
      virtual rvec CalcGridded(const double x_meas, const double y_meas,
          const double z_meas, const ThreeDGravityModel &Model,
          rmat &Sensitivities);
    public:
      TensorOMPGravityImp();
      virtual ~TensorOMPGravityImp();
      };

  }

#endif /* TENSOROMPGRAVITYIMP_H_ */
