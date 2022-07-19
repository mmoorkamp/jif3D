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
      const double threshold = 1e-30;
      //! The number of elements we expect to transform, this is purely used for potential error checking between creating the object and using it
      size_t ntrans;
      friend class access;
      //! Provide serialization to be able to store objects and, more importantly for simpler MPI parallelization
      template<class Archive>
      void serialize(Archive &ar, const unsigned int version)
        {
          ar & base_object<VectorTransform>(*this);
          ar & ntrans;
        }
    public:
      virtual size_t GetInputSize() const override
        {
          return ntrans;
        }
      virtual size_t GetOutputSize() const override
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

      virtual jif3D::rvec Transform(const jif3D::rvec &InputVector) const override
        {
          jif3D::rvec Result(2);
          double absval = std::abs(std::complex<double>(InputVector(0), InputVector(1)));
          //for impedances any magnitude less than this is negligible
          //so we use it as a minimum to avoid taking the log of zero
          if (absval < threshold)
            Result(0) = std::log(threshold);
          else
            Result(0) = std::log(absval);
          //atan2 is protected against small values, so we do not need special precautions here
          Result(1) = std::atan2(InputVector(1), InputVector(0));
          return Result;
        }

      virtual jif3D::rmat Derivative(const jif3D::rvec &InputVector) const override
        {
          double p = std::atan2(InputVector(1), InputVector(0));
          double az = std::max(threshold,
              std::abs(std::complex<double>(InputVector(0), InputVector(1))));
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
      void serialize(Archive &ar, const unsigned int version)
        {
          ar & base_object<VectorTransform>(*this);
          ar & ntrans;
        }
    public:
      virtual size_t GetInputSize() const override
        {
          return ntrans;
        }
      virtual size_t GetOutputSize() const override
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

      virtual jif3D::rvec Transform(const jif3D::rvec &InputVector) const override
        {
          jif3D::rvec Result(2);
          Result(0) = InputVector(1);
          Result(1) = -InputVector(0);
          return Result;
        }

      virtual jif3D::rmat Derivative(const jif3D::rvec &InputVector) const override
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

