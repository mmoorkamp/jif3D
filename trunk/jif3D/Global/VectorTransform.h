//============================================================================
// Name        : VectorTransform.h
// Author      : Apr 22, 2009
// Version     :
// Copyright   : 2009, mmoorkamp
//============================================================================

#ifndef VECTORTRANSFORM_H_
#define VECTORTRANSFORM_H_

#include "../Global/Serialization.h"
#include "../Global/VecMat.h"
#include "../Global/FatalException.h"
#include <cassert>

#include "Jif3DGlobal.h"

/*! \file VectorTransform.h
 * Provide function objects that transform one vector to another. The main purpose is to
 * apply transformation to data within an inversion, e.g. impedance to apparent resistivity and phase
 * or FTG tensor to an invariant. There are also transformation for models that work slightly differently
 * in the file ModelTransforms.h .
 */

namespace jif3D
  {
    /*! \addtogroup util General utility routines
     */
    /* @{ */

    //! Transform an input vector to an output vector and calculate the associated derivate
    /*! For the joint inversion we have calculate derived quantities from our data, i.e. apparent
     * resistivity for MT or an invariant for FTG data. We need an interface for these kind of transformations
     * which is provided by this base class. The transforms assume that GetInputSize elements in the input vector
     * form a logical block of data, e.g. the 9 elements of the FTG tensor and that input vector can be partioned
     * into N segments of this size. Then the output will be a vector of GetOutputSize * N elements where the
     * GetOutputSize consecutive elements correspond to one logical block in the input vector.
     *
     * For example to calculate an invariant from FTG data, GetInputSize would be 9 and GetOutputSize would give 1. An InputVector
     * of size 90, i.e. 10 tensor observations, would result in and output of 10, i.e. 10 invariants.
     */
    class J3DEXPORT VectorTransform
      {
    public:
      //! Provide serialization to be able to store objects and, more importantly for simpler MPI parallelization
      template<class Archive>
      void serialize(Archive &ar, const unsigned int version)
        {
        }
      //! How many consecutive elements in the input vector form a logical block of data that this transform works on
      virtual size_t GetInputSize() const
        {
          return 0;
        }
      //! How many elements will one logical input block be transformed to
      virtual size_t GetOutputSize() const
        {
          return 0;
        }
      //! Transform the input vector
      virtual jif3D::rvec Transform(const jif3D::rvec &InputVector) const
        {
          throw jif3D::FatalException("Illegal call to base class");
        }
      //! Give the matrix of partial derivatives with respect to the input parameters for the transformation \f$ \partial f/\partial m_i \f$.
      virtual jif3D::rmat Derivative(const jif3D::rvec &InputVector) const
        {
          throw jif3D::FatalException("Illegal call to base class");
        }
      VectorTransform()
        {
        }
      virtual ~VectorTransform()
        {
        }
      };

    //! In some cases we just want the output be the same as the input, this class simply copies the input
    /*! This class implements the simplest case, we just copy the input to the output.
     */
    class J3DEXPORT CopyTransform: public VectorTransform
      {
    private:
      //! The number of elements we expect to transform, this is purely used for potential error checking between creating the object and using it
      size_t ntrans;
    public:
      //! Provide serialization to be able to store objects and, more importantly for simpler MPI parallelization
      template<class Archive>
      void serialize(Archive &ar, const unsigned int version)
        {
          ar & ntrans;
        }
      virtual size_t GetInputSize() const override
        {
          return ntrans;
        }
      virtual size_t GetOutputSize() const override
        {
          return ntrans;
        }
      //! This transform only copies its input, but we can specify an expected size to perform error checking in code that uses the transform
      CopyTransform(size_t intendedsize) :
          ntrans(intendedsize)
        {
        }
      virtual ~CopyTransform()
        {
        }
      //! This "transformation" just passes the input parameter through
      virtual jif3D::rvec Transform(const jif3D::rvec &InputVector) const override
        {
          return InputVector;
        }
      //! When generalized and physical parameters are the same the derivative is 1 for all parameters
      virtual jif3D::rmat Derivative(const jif3D::rvec &InputVector) const override
        {
          return ublas::identity_matrix<double>(InputVector.size());
        }
      };
    //! Apply a transformation to a vector
    /*! The VectorTransform class works on a vector which length is
     * determined by the number of elements for a single transformation, i.e. 9 input
     * elements and 1 output element to calculate the invariant from the FTG tensor.
     * If we have a vector that consists of several data to which we want to apply
     * the same transformation we can use this function.
     * @param InputVector The vector containing the input data for the transformation, has to have a length of a multiple of the transformation size
     * @param Transform The transformation class to apply to the data
     * @return A vector that contains the transformed data
     */
    template<class VectorTransform> J3DEXPORT
    std::vector<double> ApplyTransform(const std::vector<double> &InputVector,
        VectorTransform &Transform)
      {
        const size_t insize = InputVector.size();
        const size_t step = Transform.GetInputSize();
        const size_t nout = Transform.GetOutputSize();
        if (insize % step != 0)
          throw jif3D::FatalException(
              "Transformation size needs to be integer multiple of input data size",
              __FILE__, __LINE__);
        std::vector<double> Output(insize / step * nout);

        jif3D::rvec curr(nout), in(step);
        for (size_t i = 0; i < insize; i += step)
          {
            const size_t start = i / step * nout;
            std::copy(InputVector.begin() + i, InputVector.begin() + i + step,
                in.begin());
            curr = Transform.Transform(in);
            std::copy(curr.begin(), curr.end(), Output.begin() + start);
          }
        return Output;
      }

    template<class VectorTransform> J3DEXPORT
    jif3D::rvec ApplyTransform(const jif3D::rvec &InputVector,
        VectorTransform &Transform)
      {
        const size_t insize = InputVector.size();
        const size_t step = Transform.GetInputSize();
        const size_t nout = Transform.GetOutputSize();
        if (insize % step != 0)
          throw jif3D::FatalException(
              "Transformation size needs to be integer multiple of input data size",
              __FILE__, __LINE__);
        jif3D::rvec Output(insize / step * nout);

        jif3D::rvec curr(nout), in(step);
        for (size_t i = 0; i < insize; i += step)
          {
            const size_t start = i / step * nout;
            std::copy(InputVector.begin() + i, InputVector.begin() + i + step,
                in.begin());
            curr = Transform.Transform(in);
            std::copy(curr.begin(), curr.end(), Output.begin() + start);
          }
        return Output;
      }

  /* @} */
  }

#endif /* VECTORTRANSFORM_H_ */
