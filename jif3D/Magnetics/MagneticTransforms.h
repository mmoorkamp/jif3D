/*
 * MagneticTransforms.h
 *
 *  Created on: Nov 19, 2013
 *      Author: mmoorkamp
 */

#ifndef MAGNETICTRANSFORMS_H_
#define MAGNETICTRANSFORMS_H_

#include "../Global/Serialization.h"

#include "../Global/VectorTransform.h"
#include "../Global/Jif3DGlobal.h"
#include <cassert>

namespace jif3D
  {
    class J3DEXPORT TotalField: public VectorTransform
      {
    private:
      //! we expect 3 magnetic field elements at each call
      static const size_t ninput = 3;
      //! we return one invariant as output at each call
      static const size_t noutput = 1;
      //! calculate the invariant from 3 observations
      double CalcTotalField(const jif3D::rvec &Data) const
        {
          return sqrt(Data(0) * Data(0) + Data(1) * Data(1) + Data(2) * Data(2));
        }
    public:
      template<class Archive>
      void serialize(Archive &ar, const unsigned int version)
        {
          ar & base_object<VectorTransform>(*this);
        }
      //! Return the size of the input vector this class expects
      virtual size_t GetInputSize() const override
        {
          return ninput;
        }
      //! Return the size of the input vector this class will yield
      virtual size_t GetOutputSize() const override
        {
          return noutput;
        }
      //! Take a vector of 9 tensor elements and calculate the invariant
      /*! This function performs the transformation of a single tensor.
       * @param InputVector The tensor elements as a vector in c-storage order, has to have 9 elements
       * @return A vector with a single element, the calculated invariant.
       */
      virtual jif3D::rvec Transform(const jif3D::rvec &InputVector) const override
        {
          assert(InputVector.size() == ninput);
          jif3D::rvec result(1);
          result(0) = CalcTotalField(InputVector);
          return result;
        }
      //! Calculate the partial derivative of the invariant with respect to the magnetic field
      /*! For a magnetic field measurement this function calculates the partial derivatives of the invariant with respect
       * to the magnetic field elements
       * @param InputVector The magnetic field elements as a vector in c-storage order, has to have 3 elements
       * @return A 1x3 matrix of partial derivatives
       */
      virtual jif3D::rmat Derivative(const jif3D::rvec &InputVector) const override
        {
          assert(InputVector.size() == ninput);
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

    class TotalFieldAnomaly: public VectorTransform
      {
    private:
      const double inclination;
      const double declination;
      double BxComp;
      double ByComp;
      double BzComp;
      //! we expect 3 magnetic field elements at each call
      static const size_t ninput = 3;
      //! we return one invariant as output at each call
      static const size_t noutput = 1;
    public:
      template<class Archive>
      void serialize(Archive &ar, const unsigned int version)
        {
          ar & base_object<VectorTransform>(*this);
          ar & inclination;
          ar & declination;
          ar & BxComp;
          ar & ByComp;
          ar & BzComp;
        }
      //! Return the size of the input vector this class expects
      virtual size_t GetInputSize() const override
        {
          return ninput;
        }
      //! Return the size of the input vector this class will yield
      virtual size_t GetOutputSize() const override
        {
          return noutput;
        }
      //! Take a vector of 3 tensor elements and calculate the total field anomaly
      /*! This function performs the transformation of a single magnetic field vector.
       * @param InputVector The 3 magnetic field components as a vector in c-storage order, has to have 3 elements
       * @return A vector with a single element, the calculated total field anomaly.
       */
      virtual jif3D::rvec Transform(const jif3D::rvec &InputVector) const override
        {
          assert(InputVector.size() == ninput);
          jif3D::rvec result(1);
          result(0) = (BxComp * InputVector(0) + ByComp * InputVector(1)
              + BzComp * InputVector(2));
          return result;
        }
      //! Calculate the partial derivative of the total field anomaly with respect to the magnetic field
      /*! For a magnetic field measurement this function calculates the partial derivatives of the total field with respect
       * to the magnetic field elements
       * @param InputVector The magnetic field elements as a vector in c-storage order, has to have 3 elements
       * @return A 1x3 matrix of partial derivatives
       */
      virtual jif3D::rmat Derivative(const jif3D::rvec &InputVector) const override
        {
          assert(InputVector.size() == ninput);
          jif3D::rmat InvarSens(noutput, ninput);

          InvarSens(0, 0) = BxComp;
          InvarSens(0, 1) = ByComp;
          InvarSens(0, 2) = BzComp;
          return InvarSens;
        }
      TotalFieldAnomaly(double inc = 0.0, double dec = 0.0) :
          inclination(inc), declination(dec)
        {
          BxComp = cos(inclination) * cos(declination);
          ByComp = cos(inclination) * sin(declination);
          BzComp = sin(inclination);
        }
      virtual ~TotalFieldAnomaly()
        {
        }
      };
  }

#endif /* MAGNETICTRANSFORMS_H_ */
