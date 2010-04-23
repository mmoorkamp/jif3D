//============================================================================
// Name        : CurvatureRegularization.h
// Author      : Jan 11, 2010
// Version     : 
// Copyright   : 2010, mmoorkamp
//============================================================================


#ifndef CURVATUREREGULARIZATION_H_
#define CURVATUREREGULARIZATION_H_

#include "MatOpRegularization.h"

namespace jiba
  {
    //! Minimize the laplacian of the Model
    /*! This objective function can be used to minimize the laplacian over a 3D model in order
     * to regularize the inversion. For each cell that is not on the model boundary of the model
     * domain we calculate \f$ \partial^2/\partial x m + \partial^2/\partial y m + \partial^2/\partial z m \f$. If
     * we specify a reference model (see description for MatOpRegularization), this reference model
     * is substracted before the calculation of the laplacian.
     */
    class CurvatureRegularization: public jiba::MatOpRegularization
      {
    private:
      //! A small number to stabilize the matrix operator
      const double Eps;
      void ConstructOperator(const jiba::ThreeDModelBase &ModelGeometry);
    public:
      //! The clone function provides a virtual constructor
      virtual CurvatureRegularization *clone() const
        {
          return new CurvatureRegularization(*this);
        }
      //! We have to provide the model geometry to the constructor
      /*! We have to specify the geometry of the 3D model by providing
       * a ThreeDModelBase or derived class. Also the gradient regularization
       * contains a small component to minimize the absolute value of the model
       * vector. This stabilizes the Operator matrix that otherwise would be singular.
       * Under normal circumstances the default value should work well, but it can be specified
       * to change the behaviour of the regularization.
       * @param Geometry A 3D model class that describes the cell geometries
       * @param ModEps The weight of the absolute value minimization of the model vector.
       */
      explicit CurvatureRegularization(const jiba::ThreeDModelBase &Geometry,
          const double ModEps = 1e-8) :
        MatOpRegularization(Geometry), Eps(ModEps)
        {
          ConstructOperator(Geometry);
        }
      virtual ~CurvatureRegularization();
      };

  }

#endif /* CURVATUREREGULARIZATION_H_ */
