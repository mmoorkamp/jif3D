//============================================================================
// Name        : GradientRegularization.h
// Author      : Jun 5, 2009
// Version     :
// Copyright   : 2009, mmoorkamp
//============================================================================


#ifndef GRADIENTREGULARIZATION_H_
#define GRADIENTREGULARIZATION_H_

#include "MatOpRegularization.h"

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
    class GradientRegularization: public jiba::MatOpRegularization
      {
    private:
      const double Eps;
      void ConstructOperator(const jiba::ThreeDModelBase &ModelGeometry);
    public:
      //! We have to provide the model geometry to the constructor
      explicit GradientRegularization(const jiba::ThreeDModelBase &Geometry, const double ModEps = 1e-8) :
        MatOpRegularization(Geometry),
        Eps(ModEps)
        {
        ConstructOperator(Geometry);
        }
      virtual ~GradientRegularization();
      };
  /* @} */
  }

#endif /* GRADIENTREGULARIZATION_H_ */
