//============================================================================
// Name        : CrossGradient.h
// Author      : Oct 22, 2009
// Version     :
// Copyright   : 2009, mmoorkamp
//============================================================================

#ifndef CROSSGRADIENT_H_
#define CROSSGRADIENT_H_

#include "../Global/Serialization.h"
#include "../Global/Jif3DGlobal.h"
#include "GradientRegularization.h"
#include "../Inversion/ObjectiveFunction.h"
#include "../ModelBase/ThreeDModelBase.h"

namespace jif3D
  {
    /** \addtogroup Regularization classes to regularize the inversion */
    /* @{ */
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
    class J3DEXPORT CrossGradient: public jif3D::ObjectiveFunction
      {
    private:
      //! The object to calculate the spatial gradient for the first model
      GradientRegularization FirstGradient;
      //! The object to calculate the spatial gradient for the second model
      GradientRegularization SecondGradient;
      //! The model geometry for the two models, this is always identical to the geometry in GradientRegularization objects
      ThreeDModelBase ModelGeometry;
      friend class access;
      //! Provide serialization to be able to store objects and, more importantly for simpler MPI parallelization
      template<class Archive>
      void serialize(Archive & ar, const unsigned int version)
        {
          ar & base_object<ObjectiveFunction>(*this);
          ar & FirstGradient;
          ar & SecondGradient;
          ar & ModelGeometry;
        }
    public:
      //! The clone function provides a virtual constructor
      virtual CrossGradient *clone() const override
        {
          return new CrossGradient(*this);
        }
      //! We never want to terminate the inversion because the regularization has reached a particular value, so we return 0
      virtual double ConvergenceLimit() const override
        {
          return -1.0;
        }
      //! The implementation of the cross-gradient calculation
      virtual void
      ImplDataDifference(const jif3D::rvec &Model, jif3D::rvec &Diff) override;
      //! The gradient of the cross-gradient objective function with respect to the model parameters
      virtual jif3D::rvec ImplGradient(const jif3D::rvec &Model, const jif3D::rvec &Diff)
          override;
      //!The constructor takes any 3D model object as a parameter to extract the geometry information
      /*! The geometry information of the 3D model is passed to the cross-gradient function through
       * the constructor. Note that during the inversion we assume Model.size() == 2*Geometry.GetData().num_elements()
       * as explained in the general description of this class.
       * @param Geometry An object that contains the information about the model geometry (cell sizes etc.)
       */
      explicit CrossGradient(const jif3D::ThreeDModelBase &Geometry) :
          FirstGradient(Geometry, 0.0), SecondGradient(Geometry, 0.0), ModelGeometry(
              Geometry)
        {

        }
      CrossGradient(const jif3D::ThreeDModelBase &Geometry,
          const jif3D::ThreeDModelBase &TearModelX,
          const jif3D::ThreeDModelBase &TearModelY,
          const jif3D::ThreeDModelBase &TearModelZ) :
          FirstGradient(Geometry, TearModelX, TearModelY, TearModelZ, 0.0), SecondGradient(
              Geometry, TearModelX, TearModelY, TearModelZ, 0.0), ModelGeometry(Geometry)
        {

        }
      virtual ~CrossGradient()
        {
        }
      };
  /* @} */
  }

#endif /* CROSSGRADIENT_H_ */
