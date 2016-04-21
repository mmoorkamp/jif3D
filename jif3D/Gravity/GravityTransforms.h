//============================================================================
// Name        : GravityTransforms.h
// Author      : Apr 22, 2009
// Version     :
// Copyright   : 2009, mmoorkamp
//============================================================================

#ifndef GRAVITYTRANSFORMS_H_
#define GRAVITYTRANSFORMS_H_

#include "../Global/VectorTransform.h"
#include "../Global/FatalException.h"

namespace jif3D
  {
    /** \addtogroup gravity Gravity forward modeling, display and inversion */
    /* @{ */

    //! Calculate I2 as described in Pedersen and Rasmussen 1990 and the associated derivative
    /*! This transformation class takes the 9 elements of the FTG tensor as input and gives
     * the invariant I2 and its partial derivatives with respect to the tensor elements as output.
     */
    class FTGInvariant: public VectorTransform
      {
    private:
      //! we expect 9 tensor elements at each call
      static const size_t ninput = 9;
      //! we return one invariant as output at each call
      static const size_t noutput = 1;
      //! calculate the invariant from 9 observations
      double CalcInvariant(const jif3D::rvec &Data) const
        {
          return Data(0) * Data(4) + Data(4) * Data(8) + Data(0) * Data(8)
              - Data(3) * Data(1) - Data(7) * Data(5) - Data(2) * Data(6);
        }
    public:
      //! Return the size of the input vector this class expects
      virtual size_t GetInputSize() const
        {
          return ninput;
        }
      //! Return the size of the input vector this class will yield
      virtual size_t GetOutputSize() const
        {
          return noutput;
        }
      //! Take a vector of 9 tensor elements and calculate the invariant
      /*! This function performs the transformation of a single tensor.
       * @param InputVector The tensor elements as a vector in c-storage order, has to have 9 elements
       * @return A vector with a single element, the calculated invariant.
       */
      virtual jif3D::rvec Transform(const jif3D::rvec &InputVector) const
        {
          assert(InputVector.size() == ninput);
          jif3D::rvec result(1);
          result(0) = CalcInvariant(InputVector);
          return result;
        }
      //! Calculate the partial derivative of the invariant with respect to the tensor elements
      /*! For a single tensor this function calculates the partial derivatives of the invariant with respect
       * to the tensor elements
       * @param InputVector The tensor elements as a vector in c-storage order, has to have 9 elements
       * @return A 1x9 matrix of partial derivatives
       */
      virtual jif3D::rmat Derivative(const jif3D::rvec &InputVector) const
        {
          const size_t ndata = InputVector.size();
          if (ndata != ninput)
            {
              throw jif3D::FatalException(
                  "Number of data does not match expected input size ! ", __FILE__, __LINE__);
            }

          jif3D::rmat InvarSens(noutput, ninput);

          InvarSens(0, 0) = InputVector(4) + InputVector(8);
          InvarSens(0, 1) = -InputVector(3);
          InvarSens(0, 2) = -InputVector(6);
          InvarSens(0, 3) = -InputVector(1);
          InvarSens(0, 4) = InputVector(0) + InputVector(8);
          InvarSens(0, 5) = -InputVector(7);
          InvarSens(0, 6) = -InputVector(2);
          InvarSens(0, 7) = -InputVector(5);
          InvarSens(0, 8) = InputVector(0) + InputVector(4);
          return InvarSens;
        }
      FTGInvariant()
        {
        }
      virtual ~FTGInvariant()
        {
        }
      };
  /* @} */
  }
#endif /* GRAVITYTRANSFORMS_H_ */
