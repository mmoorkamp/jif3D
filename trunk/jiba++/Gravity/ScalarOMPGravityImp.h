//============================================================================
// Name        : ScalarOMPGravityImp.h
// Author      : Dec 10, 2008
// Version     :
// Copyright   : 2008, mmoorkamp
//============================================================================


#ifndef SCALAROMPGRAVITYIMP_H_
#define SCALAROMPGRAVITYIMP_H_

#include "ThreeDGravityImplementation.h"

namespace jiba
  {

    class ScalarOMPGravityImp: public jiba::ThreeDGravityImplementation
      {
    private:
      static const size_t ndatapermeas = 1;
    public:
      virtual size_t GetDataPerMeasurement()
        {
          return ndatapermeas;
        }
      virtual rvec CalcBackground(const double xmeas, const double ymeas,
          const double zmeas, const double xwidth, const double ywidth,
          const double zwidth, const ThreeDGravityModel &Model,
          ublas::matrix_range<rmat> &sSensitivities);
      virtual rvec CalcGridded(const double x_meas, const double y_meas,
          const double z_meas, const ThreeDGravityModel &Model,
          ublas::matrix_range<rmat> &Sensitivities);
      ScalarOMPGravityImp();
      virtual ~ScalarOMPGravityImp();
      };

  }

#endif /* SCALAROMPGRAVITYIMP_H_ */
