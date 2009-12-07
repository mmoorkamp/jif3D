//============================================================================
// Name        : GradientRegularization.h
// Author      : Jun 5, 2009
// Version     :
// Copyright   : 2009, mmoorkamp
//============================================================================


#ifndef GRADIENTREGULARIZATION_H_
#define GRADIENTREGULARIZATION_H_

#include "ObjectiveFunction.h"
#include "../ModelBase/ThreeDModelBase.h"

namespace jiba
  {
    /** \addtogroup inversion General routines for inversion */
    /* @{ */
    //! Minimize the first derivative between the cells of a 3D Model
    /*! This class provides regularization in terms of the first derivative
     * of the model. We can provide a reference model. Then the regularization
     * will minimize the gradient of the difference between reference and current model.
     * If we don't provide a reference model it will minimize the variation of the current model.
     *
     * Note that we use a sparse matrix representation for the first derivative operator
     * which is extremely slow in debug mode when all operations are checked, but fast
     * with optimizations. Also at the moment we do not consider varying cell sizes
     * in the calculation of the model roughness.
     */
    class GradientRegularization: public jiba::ObjectiveFunction
      {
    private:
      //The operator matrix for the first spatial derivative is sparse
      //so we provide a typedef for convenience
      typedef boost::numeric::ublas::compressed_matrix<double,
          boost::numeric::ublas::column_major> comp_mat;
      //we have individual weights for each direction
      double xweight;
      double yweight;
      double zweight;
      //and store the result in a sparse matrix
      comp_mat XOperatorMatrix;
      comp_mat YOperatorMatrix;
      comp_mat ZOperatorMatrix;
      jiba::rvec XGrad;
      jiba::rvec YGrad;
      jiba::rvec ZGrad;
      // the vector of the reference model
      jiba::rvec Reference;
      //! Construct the operator matrix for the given model
      /*! From the geometry information in the ModelGeometry parameter, construct
       * the operator matrix for the first derivative
       * @param ModelGeometry The object containing the number of cells in each direction
       * @return A sparse matrix representation of the spatial first derivative
       */
      void SmoothnessOperator(const jiba::ThreeDModelBase &ModelGeometry);
      //! The misfit, i.e. roughness of the current model
      virtual void ImplDataDifference(const jiba::rvec &Model, jiba::rvec &Diff)
        {
          const size_t nmod = Model.size();
          if (Reference.size() != nmod)
            {
              Reference.resize(nmod);
              Reference.clear();
            }
          assert(XOperatorMatrix.size1() == Model.size());
          assert(YOperatorMatrix.size1() == Model.size());
          assert(ZOperatorMatrix.size1() == Model.size());
          Diff.resize(3 * nmod);
          ublas::vector_range<jiba::rvec> xrange(Diff, ublas::range(0,
              nmod));
          ublas::vector_range<jiba::rvec> yrange(Diff, ublas::range(
              nmod, 2 * nmod));
          ublas::vector_range<jiba::rvec> zrange(Diff, ublas::range(2
              * nmod, 3 * nmod));
          jiba::rvec x(Model - Reference);
          ublas::axpy_prod(XOperatorMatrix, x, xrange, true);
          ublas::axpy_prod(YOperatorMatrix, x, yrange, true);
          ublas::axpy_prod(ZOperatorMatrix, x, zrange, true);
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
          ublas::vector_range<const jiba::rvec> yrange(Diff, ublas::range(
              nmod, 2 * nmod));
          ublas::vector_range<const jiba::rvec> zrange(Diff, ublas::range(2
              * nmod, 3 * nmod));
          ublas::axpy_prod(ublas::trans(XOperatorMatrix), xrange, XGrad);
          ublas::axpy_prod(ublas::trans(YOperatorMatrix), yrange, YGrad);
          ublas::axpy_prod(ublas::trans(ZOperatorMatrix), zrange, ZGrad);
          return 2.0*(sqrt(xweight) * XGrad + sqrt(yweight) * YGrad + sqrt(zweight)*ZGrad);
        }
    public:
      void SetXWeight(const double Weight)
        {
          xweight = Weight;
        }
      void SetYWeight(const double Weight)
        {
          yweight = Weight;
        }
      void SetZWeight(const double Weight)
        {
          zweight = Weight;
        }
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
      //! We have to provide the model geometry to the constructor
      explicit GradientRegularization(const jiba::ThreeDModelBase &Geometry) :
        xweight(1.0), yweight(1.0), zweight(1.0), XOperatorMatrix(
            Geometry.GetNModelElements(), Geometry.GetNModelElements()),
            YOperatorMatrix(Geometry.GetNModelElements(),
                Geometry.GetNModelElements()), ZOperatorMatrix(
                Geometry.GetNModelElements(), Geometry.GetNModelElements())
        {
          SmoothnessOperator(Geometry);
        }
      virtual ~GradientRegularization();
      };
  /* @} */
  }

#endif /* GRADIENTREGULARIZATION_H_ */
