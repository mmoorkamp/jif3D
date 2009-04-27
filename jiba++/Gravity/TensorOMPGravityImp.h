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
    /** \addtogroup gravity Gravity forward modeling, display and inversion */
    /* @{ */
    //! Calculate a FTG gravity response using OpenMP parallelization
    /*! This class is the tensorial counterpart to ScalarOMPGravityImp.
     * It calculates the 9 elements of second derivatives of the gravitational potential.
     * It only implements the calculation of the background and the gridded part.
     * The assembly of the two parts is performed by the base class ThreeDGravityImplementation.
     */
    class TensorOMPGravityImp: public jiba::ThreeDGravityImplementation
      {
    private:
      //! Implement the calculation of the background response
      virtual rvec CalcBackground(const size_t measindex, const double xwidth, const double ywidth,
          const double zwidth, const ThreeDGravityModel &Model,
          rmat &Sensitivities);
      //! Calculate the response of the gridded part
      virtual rvec CalcGridded(const size_t measindex, const ThreeDGravityModel &Model,
          rmat &Sensitivities);
      static const size_t ndatapermeas = 9;
    public:
      virtual size_t RawDataPerMeasurement()
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
