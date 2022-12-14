/*
 * MTTransforms.h
 *
 *  Created on: 20 Apr 2016
 *      Author: mm489
 */

#include "../Global/Serialization.h"
#include "../Global/VectorTransform.h"

namespace jif3D
  {

    class ComplexLogTransform: public VectorTransform
      {
    private:
      //! The number of elements we expect to transform, this is purely used for potential error checking between creating the object and using it
      size_t ntrans;
      friend class access;
      //! Provide serialization to be able to store objects and, more importantly for simpler MPI parallelization
      template<class Archive>
      void serialize(Archive & ar, const unsigned int version)
        {
          ar & base_object<VectorTransform>(*this);
          ar & ntrans;
        }
    public:
      virtual size_t GetInputSize() const
        {
          return ntrans;
        }
      virtual size_t GetOutputSize() const
        {
          return ntrans;
        }

      ComplexLogTransform() :
          ntrans(2)
        {
        }
      virtual ~ComplexLogTransform()
        {
        }

      virtual jif3D::rvec Transform(const jif3D::rvec &InputVector) const
        {
          jif3D::rvec Result(2);
          double absval = std::abs(std::complex<double>(InputVector(0), InputVector(1)));
          Result(0) = 0.0;
          if (absval != 0.0)
            Result(0) = std::log(absval);
          Result(1) = std::atan2(InputVector(1), InputVector(0));
          return Result;
        }

      virtual jif3D::rmat Derivative(const jif3D::rvec &InputVector) const
        {
          double p = std::atan2(InputVector(1), InputVector(0));
          double az = std::abs(std::complex<double>(InputVector(0), InputVector(1)));
          jif3D::rmat Result(2, 2, 0.0);
          Result(0, 0) = cos(p) / az;
          Result(0, 1) = sin(p) / az;
          Result(1, 0) = -sin(p) / az;
          Result(1, 1) = cos(p) / az;
          return Result;
        }
      };

    class SwapTransform: public VectorTransform
      {
    private:
      //! The number of elements we expect to transform, this is purely used for potential error checking between creating the object and using it
      size_t ntrans;
      friend class access;
      //! Provide serialization to be able to store objects and, more importantly for simpler MPI parallelization
      template<class Archive>
      void serialize(Archive & ar, const unsigned int version)
        {
          ar & base_object<VectorTransform>(*this);
          ar & ntrans;
        }
    public:
      virtual size_t GetInputSize() const
        {
          return ntrans;
        }
      virtual size_t GetOutputSize() const
        {
          return ntrans;
        }

      SwapTransform() :
          ntrans(2)
        {
        }
      virtual ~SwapTransform()
        {
        }

      virtual jif3D::rvec Transform(const jif3D::rvec &InputVector) const
        {
          jif3D::rvec Result(2);
          Result(0) = InputVector(1);
          Result(1) = -InputVector(0);
          return Result;
        }

      virtual jif3D::rmat Derivative(const jif3D::rvec &InputVector) const
        {
          jif3D::rmat Result(2, 2, 0.0);
          Result(0, 0) = 0.0;
          Result(0, 1) = 1.0;
          Result(1, 0) = -1.0;
          Result(1, 1) = 0.0;
          return Result;
        }
      };
  }

