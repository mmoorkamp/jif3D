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
    class MatOpRegularization: public jiba::ObjectiveFunction
      {
    private:
      //The operator matrix for the first spatial derivative is sparse
      //so we provide a typedef for convenience
      typedef boost::numeric::ublas::compressed_matrix<double,
          boost::numeric::ublas::column_major> comp_mat;
      //we have individual weights for each direction
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
          ublas::vector_range<jiba::rvec> xrange(Diff, ublas::range(0, nmod));
          ublas::vector_range<jiba::rvec> yrange(Diff, ublas::range(nmod, 2
              * nmod));
          ublas::vector_range<jiba::rvec> zrange(Diff, ublas::range(2 * nmod, 3
              * nmod));
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
      jiba::rvec Reference;
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
