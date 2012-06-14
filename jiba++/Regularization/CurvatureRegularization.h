//============================================================================
// Name        : CurvatureRegularization.h
// Author      : Jan 11, 2010
// Version     : 
// Copyright   : 2010, mmoorkamp
//============================================================================

#ifndef CURVATUREREGULARIZATION_H_
#define CURVATUREREGULARIZATION_H_

#include <boost/serialization/serialization.hpp>
#include <boost/serialization/base_object.hpp>
#include "MatOpRegularization.h"
#include "MakeTearModel.h"

namespace jiba
  {
    /** \addtogroup Regularization classes to regularize the inversion */
    /* @{ */
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
      //! Construct the operator matrix to calculate the regularization
      /*! We construct three operator matrices, one for each spatial direction, to calculate the
       * curvature regularization term. This way we can weight each direction separately. We need
       * information about the model geometry and possible tears in the regularization in each direction.
       * @param ModelGeometry An object containing the geometry of the model we want to regularize
       * @param TearModelX Specify tearing in x-direction. The geometry has to match ModelGeometry. A value of 0 in a cell signifies that we want to exclude this cell from regularization in x-direction.
       * @param TearModelY Specify tearing in y-direction. The geometry has to match ModelGeometry. A value of 0 in a cell signifies that we want to exclude this cell from regularization in y-direction.
       * @param TearModelZ Specify tearing in z-direction. The geometry has to match ModelGeometry. A value of 0 in a cell signifies that we want to exclude this cell from regularization in z-direction.
       */
      void ConstructOperator(const jiba::ThreeDModelBase &ModelGeometry,
          const jiba::ThreeDModelBase &TearModelX,
          const jiba::ThreeDModelBase &TearModelY,
          const jiba::ThreeDModelBase &TearModelZ);
      friend class boost::serialization::access;
      //! Provide serialization to be able to store objects and, more importantly for simpler MPI parallelization
      template<class Archive>
      void serialize(Archive & ar, const unsigned int version)
        {
          ar & boost::serialization::base_object<MatOpRegularization>(*this);
          ar & Eps;

        }
    public:
      //! The clone function provides a virtual constructor
      virtual CurvatureRegularization *clone() const
        {
          return new CurvatureRegularization(*this);
        }
      //! This alternative constructor allows to introduce tears in the regularization in the three spatial directions
      /*! There are situations where we do not want to regularize in a certain part of the model. Either because we want
       * to facilitate discontinuities, or because a part of the inversion domain has been fixed and we do not want to
       * have this part influence the value of the regularization functional. In this case we can introduce tears in
       * the regularization by passing a model object for each spatial direction to this constructor that specifies
       * the regions of tear. Each tear model has to have the same geometry as the geometry object. If the "slowness" value
       * for the respective model object is zero, we introduce a tear between this cell and the cell with the higher index
       * (to the north/east/down) for all other values the regularization works as usual.
       * @param Geometry An object containing information about the geometry of the inversion domain
       * @param TearModelX Contains information about tears in the regularization in x-direction  (North)
       * @param TearModelY Contains information about tears in the regularization in y-direction  (East)
       * @param TearModelZ Contains information about tears in the regularization in z-direction  (Down)
       * @param ModEps The weight of the absolute value minimization of the model vector.
       */
      CurvatureRegularization(const jiba::ThreeDModelBase &Geometry,
          jiba::ThreeDModelBase &TearModelX, jiba::ThreeDModelBase &TearModelY,
          jiba::ThreeDModelBase &TearModelZ, const double ModEps = 1e-8) :
          MatOpRegularization(Geometry), Eps(ModEps)
        {
          //in debug mode we check that the geometry of the tear models
          //matches the geometry we specify for the inversion domain
          //we have to check all three tear models
          assert(TearModelX.GetData().shape()[0] == Geometry.GetData().shape()[0]);
          assert(TearModelX.GetData().shape()[1] == Geometry.GetData().shape()[1]);
          assert(TearModelX.GetData().shape()[2] == Geometry.GetData().shape()[2]);

          assert(TearModelY.GetData().shape()[0] == Geometry.GetData().shape()[0]);
          assert(TearModelY.GetData().shape()[1] == Geometry.GetData().shape()[1]);
          assert(TearModelY.GetData().shape()[2] == Geometry.GetData().shape()[2]);

          assert(TearModelZ.GetData().shape()[0] == Geometry.GetData().shape()[0]);
          assert(TearModelZ.GetData().shape()[1] == Geometry.GetData().shape()[1]);
          assert(TearModelZ.GetData().shape()[2] == Geometry.GetData().shape()[2]);

          ConstructOperator(Geometry, TearModelX, TearModelY, TearModelZ);
        }
      //! We have to provide the model geometry to the constructor when we do not want to have any tears
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
          jiba::ThreeDModelBase TearModel;
          //if we do not want tears in the regularization we temporarily
          //construct a dummy tear object that contains 1 everywhere
          //and therefore applies the normal regularization everywhere
          MakeTearModel(Geometry, TearModel);
          ConstructOperator(Geometry, TearModel, TearModel, TearModel);
        }
      virtual ~CurvatureRegularization();
      };
    /* @} */
  }

#endif /* CURVATUREREGULARIZATION_H_ */
