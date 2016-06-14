//============================================================================
// Name        : HOGradientRegularization.h
// Author      : May 5, 2014
// Version     :
// Copyright   : 2014, mmoorkamp
//============================================================================


#ifndef HOGRADIENTREGULARIZATION_H_
#define HOGRADIENTREGULARIZATION_H_

#include "../Global/Serialization.h"
#include "MatOpRegularization.h"
#include "MakeTearModel.h"

namespace jif3D
  {
  /** \addtogroup Regularization classes to regularize the inversion */
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
    class HOGradientRegularization: public jif3D::MatOpRegularization
      {
    private:
      //! A usually small value that controls the amount of absolute model vector minimization performed, see also the description of the constructor
      const double Eps;
      //! Construct the matrix operators that approximate the spatial gradient in the three directions
      /*! We construct three operator matrices, one for each spatial direction, to calculate the
       * gradient regularization term. This way we can weight each direction separately. We need
       * information about the model geometry and possible tears in the regularization in each direction.
       * @param ModelGeometry An object containing the geometry of the model we want to regularize
       * @param TearModelX Specify tearing in x-direction. The geometry has to match ModelGeometry.
       *        A value of 0 in a cell signifies that we want to exclude this cell from regularization in x-direction.
       * @param TearModelY Specify tearing in y-direction. The geometry has to match ModelGeometry.
       *        A value of 0 in a cell signifies that we want to exclude this cell from regularization in y-direction.
       * @param TearModelZ Specify tearing in z-direction. The geometry has to match ModelGeometry.
       *        A value of 0 in a cell signifies that we want to exclude this cell from regularization in z-direction.
       */
      void ConstructOperator(const jif3D::ThreeDModelBase &ModelGeometry,
          const jif3D::ThreeDModelBase &TearModelX,
          const jif3D::ThreeDModelBase &TearModelY,
          const jif3D::ThreeDModelBase &TearModelZ);
      friend class access;
      //! Provide serialization to be able to store objects and, more importantly for simpler MPI parallelization
      template<class Archive>
      void serialize(Archive & ar, const unsigned int version)
        {
          ar & base_object<MatOpRegularization>(*this);
          ar & Eps;

        }
    public:
      //! The clone function provides a virtual constructor
      virtual HOGradientRegularization *clone() const
        {
          return new HOGradientRegularization(*this);
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
      explicit HOGradientRegularization(const jif3D::ThreeDModelBase &Geometry,
          const double ModEps = 1e-8) :
        MatOpRegularization(Geometry), Eps(ModEps)
        {
          jif3D::ThreeDModelBase TearModel;
          //if we do not want tears in the regularization we temporarily
          //construct a dummy tear object that contains 1 everywhere
          //and therefore applies the normal regularization everywhere
          jif3D::MakeTearModel(Geometry, TearModel);
          ConstructOperator(Geometry, TearModel, TearModel, TearModel);
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
      HOGradientRegularization(const jif3D::ThreeDModelBase &Geometry,
          jif3D::ThreeDModelBase &TearModelX, jif3D::ThreeDModelBase &TearModelY,
          jif3D::ThreeDModelBase &TearModelZ, const double ModEps = 1e-8) :
        MatOpRegularization(Geometry), Eps(ModEps)
        {
          //in debug mode we check that the geometry of the tear models
          //matches the geometry we specify for the inversion domain
          assert(TearModelX.GetXCellSizes().size() == Geometry.GetXCellSizes().size());
          assert(TearModelX.GetYCellSizes().size() == Geometry.GetYCellSizes().size());
          assert(TearModelX.GetZCellSizes().size() == Geometry.GetZCellSizes().size());

          assert(TearModelY.GetXCellSizes().size() == Geometry.GetXCellSizes().size());
          assert(TearModelY.GetYCellSizes().size() == Geometry.GetYCellSizes().size());
          assert(TearModelY.GetZCellSizes().size() == Geometry.GetZCellSizes().size());

          assert(TearModelZ.GetXCellSizes().size() == Geometry.GetXCellSizes().size());
          assert(TearModelZ.GetYCellSizes().size() == Geometry.GetYCellSizes().size());
          assert(TearModelZ.GetZCellSizes().size() == Geometry.GetZCellSizes().size());

          ConstructOperator(Geometry, TearModelX, TearModelY, TearModelZ);
        }
      virtual ~HOGradientRegularization()
        {

        }
      };
  /* @} */
  }

#endif /* HOGRADIENTREGULARIZATION_H_ */
