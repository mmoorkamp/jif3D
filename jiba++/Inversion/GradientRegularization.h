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

    class GradientRegularization: public jiba::ObjectiveFunction
      {
    private:
      typedef boost::numeric::ublas::compressed_matrix<double,boost::numeric::ublas::column_major> comp_mat;
      comp_mat OperatorMatrix;
      comp_mat SmoothnessOperator();
      jiba::ThreeDModelBase ModelGeometry;
      jiba::rvec Reference;
      virtual void ImplDataDifference(const jiba::rvec &Model, jiba::rvec &Diff)
        {
          Diff.resize(Model.size());
          ublas::axpy_prod(OperatorMatrix,Model - Reference,Diff);
        }
      //! The abstract interface for the gradient calculation
      virtual jiba::rvec ImplGradient(const jiba::rvec &Model,
          const jiba::rvec &Diff)
        {
          return 2.0 * ublas::prod(ublas::trans(OperatorMatrix),Diff);
        }
    public:
      void SetReferenceModel(const jiba::rvec &Model)
        {
          Reference = Model;
        }
      GradientRegularization(const jiba::ThreeDModelBase &Geometry) :
        OperatorMatrix(), ModelGeometry(Geometry), Reference()
        {
          OperatorMatrix = SmoothnessOperator();
        }
      virtual ~GradientRegularization();
      };

  }

#endif /* GRADIENTREGULARIZATION_H_ */
