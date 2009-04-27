//============================================================================
// Name        : VectorTransform.h
// Author      : Apr 22, 2009
// Version     :
// Copyright   : 2009, mmoorkamp
//============================================================================


#ifndef VECTORTRANSFORM_H_
#define VECTORTRANSFORM_H_
#include "../Global/VecMat.h"
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
  /* @} */
  }

#endif /* VECTORTRANSFORM_H_ */
