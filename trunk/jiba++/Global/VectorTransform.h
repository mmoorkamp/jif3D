//============================================================================
// Name        : VectorTransform.h
// Author      : Apr 22, 2009
// Version     :
// Copyright   : 2009, mmoorkamp
//============================================================================


#ifndef VECTORTRANSFORM_H_
#define VECTORTRANSFORM_H_
#include "../Global/VecMat.h"
#include <cassert>

namespace jiba
  {
    /*! \addtogroup util General utility routines
     */
    /* @{ */

    //! Transform an input vector to an output vector and calculate the associated derivate
    /*! For the joint inversion we have calculate derived quantities from our data, i.e. apparent
     * resistivity for MT or an invariant for FTG data. Also, we might have to transform model
     * parameters, e.g. velocity to density. In both cases we have to transform the data, but also
     * calculate the derivative. This class provides the general interface for both cases.
     */
    class VectorTransform
      {
    public:
      virtual size_t GetInputSize() = 0;
      virtual size_t GetOutputSize() = 0;
      virtual jiba::rvec Transform(const jiba::rvec &InputVector) = 0;
      virtual jiba::rmat Derivative(const jiba::rvec &InputVector) = 0;
      VectorTransform()
        {
        }
      virtual ~VectorTransform()
        {
        }
      };

    //! For simple tests and when only inverting a single parameter the generalized and physical parameters will be identical, so we just have to pass them on
    /*! This class implements the simplest case when the generalized are the physical
     * parameters. So we just copy the input to the output.
     */
    class CopyTransform: public VectorTransform
      {
    private:
      const size_t ninput;
      const size_t noutput;
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
    template <class VectorTransform>
    jiba::rvec ApplyTransform(const jiba::rvec &InputVector, VectorTransform &Transform)
      {
        const size_t insize = InputVector.size();
        const size_t step = Transform.GetInputSize();
        const size_t nout = Transform.GetOutputSize();
        assert(insize % step  == 0);
        jiba::rvec Output(insize/step * nout );
        for (size_t i = 0; i < insize; i +=step)
          {
            jiba::rvec temp(Transform.Transform(ublas::vector_range<const jiba::rvec>(
                          InputVector, ublas::range(i, i + step))));
            copy(temp.begin(),temp.end(),Output.begin()+i/step*nout);
          }
        return Output;
      }
  /* @} */
  }

#endif /* VECTORTRANSFORM_H_ */
