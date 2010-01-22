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
    class CurvatureRegularization: public jiba::MatOpRegularization
      {
    private:
      const double Eps;
      void ConstructOperator(const jiba::ThreeDModelBase &ModelGeometry);
    public:
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
      CurvatureRegularization(const jiba::ThreeDModelBase &Geometry,
          const double ModEps = 1e-8) :
        MatOpRegularization(Geometry), Eps(ModEps)
        {
          ConstructOperator(Geometry);
        }
      virtual ~CurvatureRegularization();
      };

  }

#endif /* CURVATUREREGULARIZATION_H_ */
