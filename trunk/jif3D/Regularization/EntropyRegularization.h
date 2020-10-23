//============================================================================
// Name        : EntropyRegularization.h
// Author      : Oct 20, 2020
// Version     :
// Copyright   : 2020, mmoorkamp
//============================================================================

#ifndef ENTROPYREGULARIZATION_H_
#define ENTROPYREGULARIZATION_H_

#include "../Global/Serialization.h"
#include "../Global/Jif3DGlobal.h"
#include "../Global/VecMat.h"
#include "RegularizationFunction.h"
#include "../ModelBase/ThreeDModelBase.h"

namespace jif3D
  {
    /** \addtogroup Regularization classes to regularize the inversion */
    /* @{ */
    //! A class to regularize models based on entropy
    class J3DEXPORT EntropyRegularization: public jif3D::RegularizationFunction
      {
    private:
      jif3D::rvec CountsX;
      double xmax;
      double xmin;
      size_t nbins;
      /*!
       * @param Model The current model
       * @param Diff The difference vector, i.e. raw value of regularization, for each model cell.
       */
      virtual void ImplDataDifference(const jif3D::rvec &Model, jif3D::rvec &Diff)
          override;
      //! The gradient of the regularization with respect to the model parameters
      virtual jif3D::rvec ImplGradient(const jif3D::rvec &Model, const jif3D::rvec &Diff)
          override;
      friend class access;
      //! Provide serialization to be able to store objects and, more importantly for simpler MPI parallelization
      template<class Archive>
      void serialize(Archive & ar, const unsigned int version)
        {
          ar & base_object<RegularizationFunction>(*this);
          ar & xmax;
          ar & xmin;
          ar & nbins;
        }

    public:

      virtual EntropyRegularization *clone() const override
        {
          return new EntropyRegularization(*this);
        }
      //! We never want to terminate the inversion because the regularization has reached a particular value, so we return -1
      virtual double ConvergenceLimit() const override
        {
          return -1.0;
        }

      //! We have to specify the model geometry when constructing a regularization object
      /*! In order to understand the spatial relationship between the elements of the
       * model vector, we need the geometry of the 3D model. It should not change
       * during the lifetime of the object and is therefore a parameter for the constructor.
       * @param Geometry An object that describes the model geometry
       */
      EntropyRegularization(double xma, double xmi, size_t nb) :
          xmax(xma), xmin(xmi), nbins(nb)

        {
        }
      virtual ~EntropyRegularization()
        {

        }
      };
  /* @} */
  }

#endif /* ENTROPYREGULARIZATION_H_ */
