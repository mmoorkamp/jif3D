/*
 * MagneticTransforms.h
 *
 *  Created on: Nov 19, 2013
 *      Author: mmoorkamp
 */

#ifndef MAGNETICTRANSFORMS_H_
#define MAGNETICTRANSFORMS_H_

#include "../Global/VectorTransform.h"
#include <cassert>

namespace jif3D {
class TotalField: public VectorTransform
      {
    private:
      //! we expect 9 tensor elements at each call
      static const size_t ninput = 3;
      //! we return one invariant as output at each call
      static const size_t noutput = 1;
      //! calculate the invariant from 3 observations
      double CalcTotalField(const jif3D::rvec &Data)
        {
          return sqrt( Data(0) * Data(0) * + Data(1) * Data(1) + Data(2) * Data(2));
        }
    public:
	//! Return the size of the input vector this class expects
      virtual size_t GetInputSize()
        {
          return ninput;
        }
      //! Return the size of the input vector this class will yield
      virtual size_t GetOutputSize()
        {
          return noutput;
        }
      //! Take a vector of 9 tensor elements and calculate the invariant
      /*! This function performs the transformation of a single tensor.
       * @param InputVector The tensor elements as a vector in c-storage order, has to have 9 elements
       * @return A vector with a single element, the calculated invariant.
       */
      virtual jif3D::rvec Transform(const jif3D::rvec &InputVector)
        {
          assert(InputVector.size() == ninput);
          jif3D::rvec result(1);
          result(0) = CalcTotalField(InputVector);
          return result;
        }
      //! Calculate the partial derivative of the invariant with respect to the tensor elements
      /*! For a single tensor this function calculates the partial derivatives of the invariant with respect
       * to the tensor elements
       * @param InputVector The tensor elements as a vector in c-storage order, has to have 9 elements
       * @return A 1x9 matrix of partial derivatives
       */
      virtual jif3D::rmat Derivative(const jif3D::rvec &InputVector)
        {
          const size_t ndata = InputVector.size();
          assert(ndata == ninput);
          jif3D::rmat InvarSens(noutput, ninput);
          double Total = CalcTotalField(InputVector);

          InvarSens(0, 0) = InputVector(0) / Total;
          InvarSens(0, 1) = InputVector(1) / Total;
          InvarSens(0, 2) = InputVector(2) / Total;
          return InvarSens;
        }
      TotalField()
        {
        }
      virtual ~TotalField()
        {
        }
      };
}

#endif /* MAGNETICTRANSFORMS_H_ */
