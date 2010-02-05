//============================================================================
// Name        : MatOpRegularization.h
// Author      : Jan 11, 2010
// Version     : 
// Copyright   : 2010, mmoorkamp
//============================================================================


#ifndef MATOPREGULARIZATION_H_
#define MATOPREGULARIZATION_H_

#include "../Inversion/ObjectiveFunction.h"
#include "../Global/VecMat.h"
#include "../ModelBase/ThreeDModelBase.h"

namespace jiba
  {
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
      //The operator matrix for the first spatial derivative is sparse
      //so we provide a typedef for convenience
      typedef boost::numeric::ublas::compressed_matrix<double,
          boost::numeric::ublas::column_major> comp_mat;
      //The implementation of the data difference is independent
      //of the form of the matrix. We calculate the model roughness
      //for all three directions and combine the results depending
      //on the values of xweight, yweight and zweight
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
          assert(XOperatorMatrix.size1() == Model.size());
          assert(YOperatorMatrix.size1() == Model.size());
          assert(ZOperatorMatrix.size1() == Model.size());
          //the difference vector will contain the roughness values
          //in the three coordinate directions separately
          Diff.resize(3 * nmod);
          ublas::vector_range<jiba::rvec> xrange(Diff, ublas::range(0, nmod));
          ublas::vector_range<jiba::rvec> yrange(Diff, ublas::range(nmod, 2
              * nmod));
          ublas::vector_range<jiba::rvec> zrange(Diff, ublas::range(2 * nmod, 3
              * nmod));
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
      virtual jiba::rvec ImplGradient(const jiba::rvec &Model,
          const jiba::rvec &Diff)
        {
          const size_t nmod = Model.size();
          XGrad.resize(nmod);
          YGrad.resize(nmod);
          ZGrad.resize(nmod);
          ublas::vector_range<const jiba::rvec> xrange(Diff, ublas::range(0,
              nmod));
          ublas::vector_range<const jiba::rvec> yrange(Diff, ublas::range(nmod,
              2 * nmod));
          ublas::vector_range<const jiba::rvec> zrange(Diff, ublas::range(2
              * nmod, 3 * nmod));
          ublas::axpy_prod(ublas::trans(XOperatorMatrix), xrange, XGrad);
          ublas::axpy_prod(ublas::trans(YOperatorMatrix), yrange, YGrad);
          ublas::axpy_prod(ublas::trans(ZOperatorMatrix), zrange, ZGrad);
          return 2.0 * (sqrt(xweight) * XGrad + sqrt(yweight) * YGrad + sqrt(
              zweight) * ZGrad);
        }
    protected:
      //the weight of the regularization in x-direction
      double xweight;
      //the weight of the regularization in y-direction
      double yweight;
      //the weight of the regularization in z-direction
      double zweight;
      //Storage for the three operator matrices
      comp_mat XOperatorMatrix;
      comp_mat YOperatorMatrix;
      comp_mat ZOperatorMatrix;
      //the gradient of the x-direction part of the regularization functional
      jiba::rvec XGrad;
      //the gradient of the y-direction part of the regularization functional
      jiba::rvec YGrad;
      //the gradient of the z-direction part of the regularization functional
      jiba::rvec ZGrad;
      //the possible storage for a reference model
      jiba::rvec Reference;
    public:
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
      //! Access to the gradient of the regularization functional that works in x-direction
      const jiba::rvec &GetXGrad()
        {
          return XGrad;
        }
      const jiba::rvec &GetYGrad()
        {
          return YGrad;
        }
      const jiba::rvec &GetZGrad()
        {
          return ZGrad;
        }
      //! Set the reference model for the roughness calculation, this is optional
      void SetReferenceModel(const jiba::rvec &Model)
        {
          Reference = Model;
        }
      explicit MatOpRegularization(const jiba::ThreeDModelBase &Geometry) :
        xweight(1.0), yweight(1.0), zweight(1.0), XOperatorMatrix(
            Geometry.GetNModelElements(), Geometry.GetNModelElements()),
            YOperatorMatrix(Geometry.GetNModelElements(),
                Geometry.GetNModelElements()), ZOperatorMatrix(
                Geometry.GetNModelElements(), Geometry.GetNModelElements())
        {
        }
      virtual ~MatOpRegularization()
        {

        }
      };

  }

#endif /* MATOPREGULARIZATION_H_ */