//============================================================================
// Name        : VectorTransform.h
// Author      : Apr 22, 2009
// Version     :
// Copyright   : 2009, mmoorkamp
//============================================================================

#ifndef VECTORTRANSFORM_H_
#define VECTORTRANSFORM_H_

#include <boost/serialization/serialization.hpp>
#include <boost/serialization/base_object.hpp>
#include <cassert>
#include "../Global/VecMat.h"

/*! \file VectorTransform.h
 * Provide function objects that transform one vector to another. The main purpose is to
 * apply transformation to data within an inversion, e.g. impedance to apparent resistivity and phase
 * or FTG tensor to an invariant. There are also transformation for models that work slightly differently
 * in the file ModelTransforms.h .
 */

namespace jiba
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
    class VectorTransform
      {
    private:
      friend class boost::serialization::access;
      //! Provide serialization to be able to store objects and, more importantly for simpler MPI parallelization
      template<class Archive>
      void serialize(Archive & ar, const unsigned int version)
        {

        }
    public:
      //! How many consecutive elements in the input vector form a logical block of data that this transform works on
      virtual size_t GetInputSize() = 0;
      //! How many elements will one logical input block be transformed to
      virtual size_t GetOutputSize() = 0;
      //! Transform the input vector
      virtual jiba::rvec Transform(const jiba::rvec &InputVector) = 0;
      //! Give the matrix of partial derivatives with respect to the input parameters for the transformation \f$ \partial f/\partial m_i \f$.
      virtual jiba::rmat Derivative(const jiba::rvec &InputVector) = 0;
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
    class CopyTransform: public VectorTransform
      {
    private:
      size_t ninput;
      size_t noutput;
      friend class boost::serialization::access;
      //! Provide serialization to be able to store objects and, more importantly for simpler MPI parallelization
      template<class Archive>
      void serialize(Archive & ar, const unsigned int version)
        {
          ar & boost::serialization::base_object<VectorTransform>(*this);
          ar & ninput;
          ar & noutput;
        }

    public:
      virtual size_t GetInputSize()
        {
          return ninput;
        }
      virtual size_t GetOutputSize()
        {
          return noutput;
        }
      //! This transform only copies its input, but we can specify an expected size to perform error checking in code that uses the transform
      CopyTransform(size_t intendedsize = 1) :
          ninput(intendedsize), noutput(intendedsize)
        {
        }
      virtual ~CopyTransform()
        {
        }
      //! This "transformation" just passes the input parameter through
      virtual jiba::rvec Transform(const jiba::rvec &InputVector)
        {
          return InputVector;
        }
      //! When generalized and physical parameters are the same the derivative is 1 for all parameters
      virtual jiba::rmat Derivative(const jiba::rvec &InputVector)
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
    template<class VectorTransform>
    jiba::rvec ApplyTransform(const jiba::rvec &InputVector, VectorTransform &Transform)
      {
        const size_t insize = InputVector.size();
        const size_t step = Transform.GetInputSize();
        const size_t nout = Transform.GetOutputSize();
        assert(insize % step == 0);
        jiba::rvec Output(insize / step * nout);
        for (size_t i = 0; i < insize; i += step)
          {
            jiba::rvec temp(
                Transform.Transform(
                    ublas::vector_range<const jiba::rvec>(InputVector,
                        ublas::range(i, i + step))));
            copy(temp.begin(), temp.end(), Output.begin() + i / step * nout);
          }
        return Output;
      }
  /* @} */
  }

#endif /* VECTORTRANSFORM_H_ */
