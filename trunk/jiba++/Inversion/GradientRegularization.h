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
      //and store the result in a sparse matrix
      comp_mat OperatorMatrix;
      // the vector of the reference model
      jiba::rvec Reference;
      //! Construct the operator matrix for the given model
      /*! From the geometry information in the ModelGeometry parameter, construct
       * the operator matrix for the first derivative
       * @param ModelGeometry The object containing the number of cells in each direction
       * @return A sparse matrix representation of the spatial first derivative
       */
      comp_mat SmoothnessOperator(const jiba::ThreeDModelBase &ModelGeometry);
      //! The misfit, i.e. roughness of the current model
      virtual void ImplDataDifference(const jiba::rvec &Model, jiba::rvec &Diff)
        {
          assert(OperatorMatrix.size1() == Model.size());
          Diff.resize(Model.size());
          if (Reference.size() == Model.size())
            {
              ublas::axpy_prod(OperatorMatrix, Model - Reference, Diff);
            }
          else
            {
              ublas::axpy_prod(OperatorMatrix, Model, Diff);
            }
        }
      //! The gradient of the regularization with respect to the model parameters
      virtual jiba::rvec ImplGradient(const jiba::rvec &Model,
          const jiba::rvec &Diff)
        {
          return 2.0 * ublas::prod(ublas::trans(OperatorMatrix), Diff);
        }
    public:
      //! Set the reference model for the roughness calculation, this is optional
      void SetReferenceModel(const jiba::rvec &Model)
        {
          Reference = Model;
        }
      //! We have to provide the model geometry to the constructor
      explicit GradientRegularization(const jiba::ThreeDModelBase &Geometry) :
        OperatorMatrix(SmoothnessOperator(Geometry)), Reference()
        {
        }
      virtual ~GradientRegularization();
      };
  /* @} */
  }

#endif /* GRADIENTREGULARIZATION_H_ */
