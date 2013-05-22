/*
 * DensityTransform.h
 *
 *  Created on: 13.02.2013
 *      Author:  mm489
 */

#ifndef DENSITYTRANSFORM_H_
#define DENSITYTRANSFORM_H_

#include <boost/serialization/serialization.hpp>
#include <boost/serialization/base_object.hpp>
#include <boost/shared_ptr.hpp>
#include "../ModelBase/ThreeDModelBase.h"
#include "GeneralModelTransform.h"

namespace jif3D
  {
    /** \addtogroup inversion General routines for inversion */
    /* @{ */

    //! Transform generalized model parameters for slowness to density
    /*! This transformation class can be used to transform generalized
     * model parameters for slowness to density using a linear
     * relationship between velocity and density
     * \f$ \rho = (1/s + b)/a \f$. This type of transformation
     * is motivated by the common velocity-density relationships
     * used in inversion problems and also observed for the Faeroe data.
     */
    class DensityTransform: public jif3D::GeneralModelTransform
      {
    private:
      //! A pointer to a transformation that gives slowness
      boost::shared_ptr<GeneralModelTransform> SlownessTransform;
      //! An object indicating where to apply the parameter relationship (value 1)
      jif3D::ThreeDModelBase RelModel;
      //! The value to use for density where the relationship is not valid
      double replacevalue;
      //! The slope for the linear relationship
      double a;
      //! b/a is the abscissa of the linear relationship
      double b;
      friend class boost::serialization::access;
      //! Provide serialization to be able to store objects and, more importantly for simpler MPI parallelization
      template<class Archive>
      void serialize(Archive & ar, const unsigned int version)
        {
          ar & boost::serialization::base_object<GeneralModelTransform>(*this);
          ar & SlownessTransform;
          ar & RelModel;
          ar & replacevalue;
          ar & a;
          ar & b;
        }
      DensityTransform() :
          replacevalue(0), a(0), b(0)
        {
        }
    public:
      //! We setup a clone function to have a virtual constructor and create polymorphic copies
      virtual DensityTransform* clone() const
        {
          return new DensityTransform(*this);
        }
      //! Transform the normalized model parameters back to physical parameters, in this case from Slowness to Density
      virtual jif3D::rvec GeneralizedToPhysical(const jif3D::rvec &FullModel) const
        {
          assert(RelModel.GetData().num_elements() == FullModel.size());
          jif3D::rvec Slowness(SlownessTransform->GeneralizedToPhysical(FullModel));
          jif3D::rvec Density(FullModel.size());
          for (size_t i = 0; i < FullModel.size(); ++i)
            {
              //check if we should apply the transform for this grid cell
              if (RelModel.GetData().data()[i])
                {
                  Density(i) = (1.0 / Slowness(i) + b) / a;
                }
              else
                {
                  Density(i) = replacevalue;
                }

            }
          return Density;
        }
      //! Transform from Density to Slowness
      virtual jif3D::rvec PhysicalToGeneralized(const jif3D::rvec &FullModel) const
        {
          jif3D::rvec Slowness(FullModel.size());
          for (size_t i = 0; i < FullModel.size(); ++i)
            {
              //check if we should apply the transform for this grid cell
              if (RelModel.GetData().data()[i])
                {
                  Slowness(i) = 1.0 / (a * FullModel(i) - b);
                }
              else
                {
                  Slowness(i) = replacevalue;
                }
            }
          return SlownessTransform->PhysicalToGeneralized(Slowness);
        }
      //! Transform the derivative with respect to the Slowness to Density
      virtual jif3D::rvec Derivative(const jif3D::rvec &FullModel,
          const jif3D::rvec &Derivative) const
        {
          jif3D::rvec Slowness(SlownessTransform->GeneralizedToPhysical(FullModel));
          jif3D::rvec SlowDeriv(SlownessTransform->Derivative(FullModel, Derivative));
          jif3D::rvec Output(FullModel.size());
          for (size_t i = 0; i < FullModel.size(); ++i)
            {
              //check if we should apply the transform for this grid cell
              if (RelModel.GetData().data()[i])
                {
                  Output(i) = -1.0 / (Slowness(i) * Slowness(i)) / a * SlowDeriv(i);
                }
              else
                {
                  Output(i) = 0.0;
                }
            }
          return Output;
        }
      //! The constructor needs a pointer to an object that gives slowness
      /*! To reduce the amount of code in the main program this class
       * takes a pointer to a transformation object that it uses to
       * transform the generalized model parameters to slowness before
       * then transforming to density.
       * We assume a functional relationship of the form \f$ \rho = (1/s +b)/a \f$,
       * the coefficients are specified in the constructor.
       * @param SlowTrans A pointer to an object that gives slowness
       * @param RModel A model object indicating where the relationship should not be applied (value 0) and where it should be applied (value 1)
       * @param rvalue Value for density for model cells where the parameter relationship does not apply
       * @param aval The slope a of the functional relationship
       * @param bval The offset value (see equation above)
       */
      DensityTransform(boost::shared_ptr<GeneralModelTransform> SlowTrans,
          const jif3D::ThreeDModelBase &RModel, double rvalue = 0.0, double aval = 5000,
          double bval = 8500) :
          SlownessTransform(SlowTrans), RelModel(RModel), replacevalue(rvalue), a(aval), b(
              bval)
        {
        }
      virtual ~DensityTransform()
        {
        }
      };
  /* @} */
  }

#endif /* DENSITYTRANSFORM_H_ */
