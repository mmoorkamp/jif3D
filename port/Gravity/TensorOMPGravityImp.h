//============================================================================
// Name        : TensorOMPGravityImp.h
// Author      : Dec 10, 2008
// Version     :
// Copyright   : 2008, mmoorkamp
//============================================================================


#ifndef TENSOROMPGRAVITYIMP_H_
#define TENSOROMPGRAVITYIMP_H_

#include "../Global/Serialization.h"
#include "../Global/Jif3DGlobal.h"
#include "../Gravity/ThreeDGravityModel.h"
#include "../GravMag/ThreeDGravMagImplementation.h"

namespace jif3D
  {
    /** \addtogroup gravity Gravity forward modeling, display and inversion */
    /* @{ */
    //! Calculate a FTG gravity response using OpenMP parallelization
    /*! This class is the tensorial counterpart to ScalarOMPGravityImp.
     * It calculates the 9 elements of second derivatives of the gravitational potential.
     * It only implements the calculation of the background and the gridded part.
     * The assembly of the two parts is performed by the base class ThreeDGravityImplementation.
     */
    class J3DEXPORT TensorOMPGravityImp: public jif3D::ThreeDGravMagImplementation<ThreeDGravityModel>
      {
    private:
      //! Implement the calculation of the background response
      virtual rvec CalcBackground(const size_t measindex, const double xwidth,
          const double ywidth, const double zwidth,
          const ThreeDGravityModel &Model, rmat &Sensitivities);
      //! Calculate the response of the gridded part
      virtual rvec CalcGridded(const size_t measindex,
          const ThreeDGravityModel &Model, rmat &Sensitivities);
      //! The gravity tensor has 9 elements and we return all of them even though they are not all independent
      static const size_t ndatapermeas = 9;
      friend class access;
      //! Provide serialization to be able to store objects and, more importantly for simpler MPI parallelization
      template<class Archive>
      void serialize(Archive & ar, const unsigned int version)
        {
          ar & base_object<ThreeDGravMagImplementation>(*this);
        }
    public:
      //! How many data do we return before any transformation
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
