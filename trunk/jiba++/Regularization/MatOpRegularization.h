//============================================================================
// Name        : MatOpRegularization.h
// Author      : Jan 11, 2010
// Version     : 
// Copyright   : 2010, mmoorkamp
//============================================================================

#ifndef MATOPREGULARIZATION_H_
#define MATOPREGULARIZATION_H_

#include <boost/serialization/serialization.hpp>
#include <boost/serialization/base_object.hpp>
#include "../Inversion/ObjectiveFunction.h"
#include "../Global/VecMat.h"
#include "../ModelBase/ThreeDModelBase.h"

namespace jiba
  {
    /** \addtogroup Regularization classes to regularize the inversion */
    /* @{ */
    //! A base class for regularization Functionals that can be expressed as Matrix operators
    /*! The two standard forms of regularization, model gradient and laplacian of the model,
     * can be easily expressed as sparse matrix operations on the model vector. Independent
     * on the form of that matrix, the calculation of the objective function and the gradient
     * calculation is identical. This class therefore implements this functionality and derived
     * classes only have to make sure to initialize the three matrices XOperatorMatrix, YOperatorMatrix
     * and ZOperatorMatrix. These calculate the model roughness in x-direction, y-direction and z-direction, respectively.
     *
     */
    class MatOpRegularization: public jiba::ObjectiveFunction
      {
    private:
      /*! The implementation of the data difference is independent
       * of the form of the matrix. We calculate the model roughness
       * for all three directions and combine the results depending
       * on the values of xweight, yweight and zweight.
       * @param Model The current model
       * @param Diff The difference vector, i.e. raw value of regularization, for each model cell.
       */
      virtual void ImplDataDifference(const jiba::rvec &Model, jiba::rvec &Diff)
        {
          const size_t nmod = Model.size();
          //if we didn't specify a reference model that should be substracted
          //before roughness calculation, we set the vector to the right size
          //and fill it with zeros. That way we do not have to consider this special
          //case any more
          if (Reference.size() != nmod)
            {
              Reference.resize(nmod);
              Reference.clear();
            }
          //make sure that the operator matrix has the right dimension
          assert(XOperatorMatrix.size1() == nmod);
          assert(YOperatorMatrix.size1() == nmod);
          assert(ZOperatorMatrix.size1() == nmod);
          //the difference vector will contain the roughness values
          //in the three coordinate directions separately
          Diff.resize(3 * nmod);
          ublas::vector_range<jiba::rvec> xrange(Diff, ublas::range(0, nmod));
          ublas::vector_range<jiba::rvec> yrange(Diff, ublas::range(nmod, 2 * nmod));
          ublas::vector_range<jiba::rvec> zrange(Diff, ublas::range(2 * nmod, 3 * nmod));
          jiba::rvec x(Model - Reference);
          //calculate the action of the regularization operator
          //on the model vector for the three coordinate directions
          //and store the difference in the corresponding range of
          //the misfit vector
          ublas::axpy_prod(XOperatorMatrix, x, xrange, true);
          ublas::axpy_prod(YOperatorMatrix, x, yrange, true);
          ublas::axpy_prod(ZOperatorMatrix, x, zrange, true);

          //weight each section of the misfit vector appropriately
          xrange *= sqrt(xweight);
          yrange *= sqrt(yweight);
          zrange *= sqrt(zweight);
        }
      //! The gradient of the regularization with respect to the model parameters
      virtual jiba::rvec ImplGradient(const jiba::rvec &Model, const jiba::rvec &Diff)
        {
          //set the gradient vector for each individual direction to the correct size
          const size_t nmod = Model.size();
          XGrad.resize(nmod);
          YGrad.resize(nmod);
          ZGrad.resize(nmod);
          //define some ranges on the difference vector (length 3*nmod)
          //that correspond to the gradients for each spatial direction
          ublas::vector_range<const jiba::rvec> xrange(Diff, ublas::range(0, nmod));
          ublas::vector_range<const jiba::rvec> yrange(Diff,
              ublas::range(nmod, 2 * nmod));
          ublas::vector_range<const jiba::rvec> zrange(Diff,
              ublas::range(2 * nmod, 3 * nmod));
          //The gradient is simply the transpose of the operator matrix
          //times the misfit, we calculate this for each direction
          ublas::axpy_prod(ublas::trans(XOperatorMatrix), xrange, XGrad);
          ublas::axpy_prod(ublas::trans(YOperatorMatrix), yrange, YGrad);
          ublas::axpy_prod(ublas::trans(ZOperatorMatrix), zrange, ZGrad);
          //the total gradient is the weighted sum of all individual gradients
          return 2.0
              * (sqrt(xweight) * XGrad + sqrt(yweight) * YGrad + sqrt(zweight) * ZGrad);
        }
      friend class boost::serialization::access;
      //! Provide serialization to be able to store objects and, more importantly for simpler MPI parallelization
      template<class Archive>
      void serialize(Archive & ar, const unsigned int version)
        {
          ar & boost::serialization::base_object<ObjectiveFunction>(*this);
          ar & xweight;
          ar & yweight;
          ar & zweight;
          ar & XOperatorMatrix;
          ar & YOperatorMatrix;
          ar & ZOperatorMatrix;
          ar & XGrad;
          ar & YGrad;
          ar & ZGrad;
          ar & Reference;
        }
    protected:
      //! The weight of the regularization in x-direction
      double xweight;
      //! The weight of the regularization in y-direction
      double yweight;
      //! The weight of the regularization in z-direction
      double zweight;
      //! Storage for the operator matrix in x-direction
      comp_mat XOperatorMatrix;
      //! Storage for the operator matrix in y-direction
      comp_mat YOperatorMatrix;
      //! Storage for the operator matrix in z-direction
      comp_mat ZOperatorMatrix;
      //! The gradient of the x-direction part of the regularization functional
      jiba::rvec XGrad;
      //! The gradient of the y-direction part of the regularization functional
      jiba::rvec YGrad;
      //! The gradient of the z-direction part of the regularization functional
      jiba::rvec ZGrad;
      //! The storage for a reference model
      jiba::rvec Reference;
    public:
      //! We need a virtual constructor to create a new object from a pointer to a base class;
      /*! There are situations where we only have a pointer to the base class, but we need
       * a copy of the derived class without knowing what derived type it has. This virtual
       * constructor definition allows us to do this. Each derived class has to define this
       * function a return a pointer to a copy of itself.
       * @return A pointer to copy of the derived object
       */
      virtual MatOpRegularization *clone() const = 0;
      //! We never want to terminate the inversion because the regularization has reached a particular value, so we return -1
      virtual double ConvergenceLimit() const
        {
          return -1.0;
        }
      //! Set the regularization strength in x-direction, will be 1 otherwise
      void SetXWeight(const double Weight)
        {
          xweight = Weight;
        }
      //! Set the regularization strength in y-direction, will be 1 otherwise
      void SetYWeight(const double Weight)
        {
          yweight = Weight;
        }
      //! Set the regularization strength in z-direction, will be 1 otherwise
      void SetZWeight(const double Weight)
        {
          zweight = Weight;
        }
      //! Read only access to the gradient of the regularization functional that works in x-direction
      const jiba::rvec &GetXGrad() const
        {
          return XGrad;
        }
      const comp_mat &GetXOperator() const
        {
    	  return XOperatorMatrix;
        }
      const comp_mat &GetYOperator() const
        {
          return YOperatorMatrix;
        }
      const comp_mat &GetZOperator() const
        {
          return ZOperatorMatrix;
         }
      //! Read only access to the gradient of the regularization functional that works in y-direction
      const jiba::rvec &GetYGrad() const
        {
          return YGrad;
        }
      //! Read only access to the gradient of the regularization functional that works in z-direction
      const jiba::rvec &GetZGrad() const
        {
          return ZGrad;
        }
      //! Set the reference model for the roughness calculation, this is optional
      void SetReferenceModel(const jiba::rvec &Model)
        {
          Reference = Model;
        }
      //! We have to specify the model geometry when constructing a regularization object
      /*! In order to understand the spatial relationship between the elements of the
       * model vector, we need the geometry of the 3D model. It should not change
       * during the lifetime of the object and is therefore a parameter for the constructor.
       * @param Geometry An object that describes the model geometry
       */
      explicit MatOpRegularization(const jiba::ThreeDModelBase &Geometry) :
          xweight(1.0), yweight(1.0), zweight(1.0), XOperatorMatrix(
              Geometry.GetNModelElements(), Geometry.GetNModelElements()), YOperatorMatrix(
              Geometry.GetNModelElements(), Geometry.GetNModelElements()), ZOperatorMatrix(
              Geometry.GetNModelElements(), Geometry.GetNModelElements())
        {
        }
      virtual ~MatOpRegularization()
        {

        }
      };
  /* @} */
  }

#endif /* MATOPREGULARIZATION_H_ */
