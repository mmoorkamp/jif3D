//============================================================================
// Name        : CrossGradient.h
// Author      : Oct 22, 2009
// Version     : 
// Copyright   : 2009, mmoorkamp
//============================================================================


#ifndef CROSSGRADIENT_H_
#define CROSSGRADIENT_H_

#include "ObjectiveFunction.h"
#include "GradientRegularization.h"
#include "../ModelBase/ThreeDModelBase.h"
#include "../Gravity/ThreeDGravityModel.h"

namespace jiba
  {
    //! This class provides the cross-gradient constraint objective function as suggested by Gallardo and Mejo
    /*! One way to achieve a structural coupling of different models in joint inversion is through the cross-gradient constraint.
     * Given two models \f$ m_1\f$ and \f$ m_2\f$ defined on the same grid, but with different physical parameters, the cross gradient
     * is defined as \f$ \mathbf{t} = \nabla m_1 \times \nabla m_2 \f$. This function is minimal when the spatial gradients of
     * both models point in the same direction or one of them vanishes. By minimizing an objective function based on the cross-gradient
     * in the joint inversion, we enforce models with similar structural apperance.
     *
     * This class takes the grid information as a parameter to the constructor, similar to the class GradientRegularization. During the inversion
     * the model vector is interpreted as a concatenation of two models of length 2*m. The first m values are interpreted as the model vector
     * \f$ m_1\f$ that can be mapped onto the grid described by ModelGeometry. The last m values are interpreted as the corresponding values
     * in  \f$ m_2\f$.
     *
     */
    class CrossGradient: public jiba::ObjectiveFunction
      {
    private:
      GradientRegularization FirstGradient;
      GradientRegularization SecondGradient;
      ThreeDModelBase ModelGeometry;
    public:
      //! The implementation of the cross-gradient calculation
      virtual void
          ImplDataDifference(const jiba::rvec &Model, jiba::rvec &Diff);
      //! The gradient of the cross-gradient objective function with respect to the model parameters
      virtual jiba::rvec ImplGradient(const jiba::rvec &Model,
          const jiba::rvec &Diff);
      //!The constructor takes any 3D model object as a parameter to extract the geometry information
      /*! The geometry information of the 3D model is passed to the cross-gradient function through
       * the constructor. Note that during the inversion we assume Model.size() == 2*Geometry.GetData().num_elements()
       * as explained in the general description of this class.
       * @param Geometry An object that contains the information about the model geometry (cell sizes etc.)
       */
      explicit CrossGradient(const jiba::ThreeDModelBase &Geometry) :
        FirstGradient(Geometry), SecondGradient(Geometry), ModelGeometry(
            Geometry)
        {

        }
      virtual ~CrossGradient()
        {
        }
      };

  }

#endif /* CROSSGRADIENT_H_ */
