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
    /** \addtogroup gravity Gravity forward modelling, display and inversion */
    /* @{ */
    class TensorOMPGravityImp: public jiba::ThreeDGravityImplementation
      {
    private:
      virtual rvec CalcBackground(const size_t measindex, const double xwidth, const double ywidth,
          const double zwidth, const ThreeDGravityModel &Model,
          rmat &Sensitivities);
      virtual rvec CalcGridded(const size_t measindex, const ThreeDGravityModel &Model,
          rmat &Sensitivities);
      static const size_t ndatapermeas = 9;
    public:
      virtual size_t GetDataPerMeasurement()
        {
          return ndatapermeas;
        }

    public:
      TensorOMPGravityImp();
      virtual ~TensorOMPGravityImp();
      };
  /* @} */
  }

#endif /* TENSOROMPGRAVITYIMP_H_ */
