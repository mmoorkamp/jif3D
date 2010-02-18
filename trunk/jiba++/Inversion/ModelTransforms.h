//============================================================================
// Name        : ModelTransforms.h
// Author      : Apr 15, 2009
// Version     :
// Copyright   : 2009, mmoorkamp
//============================================================================

#ifndef MODELTRANSFORMS_H_
#define MODELTRANSFORMS_H_

#include "ModelDistributor.h"

/*! \file ModelTransforms.h
 * This file contains various classes to transform model parameters within an inversion.
 * These are used to either make the inversion more well behaved or to calculate
 * one physical quantity from another, like density from slowness in LogDensityTransform.
 */

namespace jiba
  {
    /** \addtogroup inversion General routines for inversion */
    /* @{ */
    //! Normalize the model parameters by dividing by a reference model.
    /*! This class takes a reference model, e.g. the starting model in the
     * inversion and divides each model parameter by the corresponding value
     * in the reference model. The two models therefore have to have the same length.
     * This makes the model parameters dimensionless and, at least in the beginning, on the
     * order of unity therefore helping to avoid problems with greatly varying magnitudes.
     */
    class NormalizeTransform: public jiba::GeneralModelTransform
      {
    private:
      // the Reference model
      const jiba::rvec Reference;
    public:
      //! Transform the normalized model parameters back to physical parameters
      virtual jiba::rvec GeneralizedToPhysical(const jiba::rvec &FullModel) const
        {
          assert(FullModel.size() == Reference.size());
          return ublas::element_prod(FullModel, Reference);
        }
      //! Transform the physical model parameters to generalized model parameters
      virtual jiba::rvec PhysicalToGeneralized(const jiba::rvec &FullModel) const
        {
          assert(FullModel.size() == Reference.size());
          return ublas::element_div(FullModel, Reference);
        }
      //! Transform the derivative with respect to the physical parameters to normalized parameters
      virtual jiba::rvec Derivative(const jiba::rvec &FullModel,
          const jiba::rvec &Derivative) const
        {
          assert(FullModel.size() == Reference.size());
          assert(FullModel.size() == Derivative.size());
          return ublas::element_prod(Derivative, Reference);
        }
      //! The constructor needs the reference model, this has to have the same size as the inversion model
      NormalizeTransform(const jiba::rvec &Ref) :
        Reference(Ref)
        {
        }
      virtual ~NormalizeTransform()
        {
        }
      };
    //! Transform normalized logarithmic parameters
    /*! This transform is used for model parameters of the form \f$m^{\star} = \ln m/m_0\f$.
     * Here \f$m_0\f$ is the reference mdoel, e.g. the starting model of the inversion.
     * The advantage over simple normalization is that we enforce the physical parameters
     * to be positive and that we can capture large variations in physical parameters
     * in a relatively small range of inversion parameters.
     */
    class LogTransform: public jiba::GeneralModelTransform
      {
    private:
      const jiba::rvec Reference;
    public:
      //! Transform the normalized model parameters back to physical parameters
      virtual jiba::rvec GeneralizedToPhysical(const jiba::rvec &FullModel) const
        {
          assert(FullModel.size() == Reference.size());
          jiba::rvec Output(FullModel.size());
          for (size_t i = 0; i < FullModel.size(); ++i)
            Output(i) = exp(FullModel(i)) * Reference(i);
          return Output;
        }
      //! Transform the physical model parameters to generalized model parameters
      virtual jiba::rvec PhysicalToGeneralized(const jiba::rvec &FullModel) const
        {
          assert(FullModel.size() == Reference.size());
          jiba::rvec Output(FullModel.size());
          for (size_t i = 0; i < FullModel.size(); ++i)
            Output(i) = log(FullModel(i) / Reference(i));
          return Output;
        }
      //! Transform the derivative with respect to the physical parameters to normalized parameters
      virtual jiba::rvec Derivative(const jiba::rvec &FullModel,
          const jiba::rvec &Derivative) const
        {

          jiba::rvec Output(FullModel.size());
          for (size_t i = 0; i < FullModel.size(); ++i)
            {
              Output(i) = Reference(i) * exp(FullModel(i)) * Derivative(i);
            }
          return Output;
        }
      //! The constructor takes a reference model for normalization as argument
      /*! When we construct the transformation object we need the values for \f$m_0\f$ as described in
       * the general description for this class.
       * @param Ref The reference model vector \f$m_0\f$
       */
      LogTransform(const jiba::rvec &Ref) :
        Reference(Ref)
        {
        }
      virtual ~LogTransform()
        {
        }
      };

    //! Change an unconstrained optimization problem to a constrained optimization problem throug a tanh transformation
    /*! We often want to constrain the range of possible values in our optimization problem between an upper and
     * a lower limit. Instead of using a constrained optimization method, we use this transformation
     * \f$ m^{\star} =  \atanh \left(2.0 * \frac{m - m_{min}} {m_{max} - m_{min}} \right) \f$ between the generalized
     * model parameters \f$ m \f$ and the generalized model parameters \f$ m^{\star} \f$. The generalized model
     * parameters can then vary over the whole range, while the physical model parameters will always been between
     * the two bounds.
     */
    class TanhTransform: public jiba::GeneralModelTransform
      {
    private:
      //! The minimum value for the pysical model parameters
      const double min;
      //! The maximum value for the pysical model parameters
      const double max;
    public:
      //! Transform the normalized model parameters back to physical parameters
      virtual jiba::rvec GeneralizedToPhysical(const jiba::rvec &FullModel) const
        {
          jiba::rvec Output(FullModel.size());
          for (size_t i = 0; i < FullModel.size(); ++i)
            Output(i) = min + (1.0 + tanh(FullModel(i))) / 2.0 * (max - min);
          return Output;
        }
      //! Transform the physical model parameters to generalized model parameters
      virtual jiba::rvec PhysicalToGeneralized(const jiba::rvec &FullModel) const
        {
          jiba::rvec Output(FullModel.size());
          for (size_t i = 0; i < FullModel.size(); ++i)
            {
              const double argument = 2.0 * (FullModel(i) - min) / (max - min)
                  - 1;
              Output(i) = atanh(argument);
            }
          return Output;
        }
      //! Transform the derivative with respect to the physical parameters to normalized parameters
      virtual jiba::rvec Derivative(const jiba::rvec &FullModel,
          const jiba::rvec &Derivative) const
        {

          jiba::rvec Output(FullModel.size());
          for (size_t i = 0; i < FullModel.size(); ++i)
            {
              Output(i) = (max - min) / (2.0 * pow(cosh(FullModel(i)), 2))
                  * Derivative(i);
            }
          return Output;
        }
      //! The constructor need the minimum and maximum physical model parameter
      /*! When we construct the transformation object we need to specify the
       * upper and lower bound for the physical model parameters.
       * @param minval The minimum value for the physical model parameters.
       * @param maxval The maximum value for the physical model parameters.
       */
      TanhTransform(const double minval = 1.0, const double maxval = 5.0) :
        min(minval), max(maxval)
        {
        }
      virtual ~TanhTransform()
        {
        }
      };

    //! Transform generalized model parameters for slowness to density
    /*! This transformation class can be used to transform generalized
     * model parameters for slowness to density using a linear
     * relationship between velocity and density. This type of transformation
     * is motivated by the common velocity-density relationships
     * used in inversion problems and also observed for the Faeroe data.
     */
    class DensityTransform: public jiba::GeneralModelTransform
      {
    private:
      //! A pointer to a transformation that gives slowness
      boost::shared_ptr<GeneralModelTransform> SlownessTransform;
    public:
      //! Transform the normalized model parameters back to physical parameters
      virtual jiba::rvec GeneralizedToPhysical(const jiba::rvec &FullModel) const
        {
          jiba::rvec Slowness(SlownessTransform->GeneralizedToPhysical(
              FullModel));
          jiba::rvec Output(FullModel.size());
          for (size_t i = 0; i < FullModel.size(); ++i)
            {
              Output(i) = (1.0 / Slowness(i) + 8500.0) / 5000.0;
            }
          return Output;
        }
      //! Transform the physical model parameters to generalized model parameters
      virtual jiba::rvec PhysicalToGeneralized(const jiba::rvec &FullModel) const
        {

          jiba::rvec Output(FullModel.size());
          for (size_t i = 0; i < FullModel.size(); ++i)
            {
              Output(i) = 1.0 / (5000.0 * FullModel(i) - 8500.0);
            }
          return SlownessTransform->PhysicalToGeneralized(Output);
        }
      //! Transform the derivative with respect to the physical parameters to normalized parameters
      virtual jiba::rvec Derivative(const jiba::rvec &FullModel,
          const jiba::rvec &Derivative) const
        {
          jiba::rvec Slowness(SlownessTransform->GeneralizedToPhysical(
              FullModel));
          jiba::rvec SlowDeriv(SlownessTransform->Derivative(FullModel,
              Derivative));
          jiba::rvec Output(FullModel.size());
          for (size_t i = 0; i < FullModel.size(); ++i)
            {
              Output(i) = -1.0 / (Slowness(i) * Slowness(i)) / 5000.0
                  * SlowDeriv(i);
            }
          return Output;
        }
      //! The constructor needs a pointer to an object that gives slowness
      /*! To reduce the amount of code in the main program this class
       * takes a pointer to a transformation object that it uses to
       * transform the generalized model parameters to slowness before
       * then transforming to density.
       * @param SlowTrans A pointer to an object that gives slowness
       */
      DensityTransform(boost::shared_ptr<GeneralModelTransform> SlowTrans) :
        SlownessTransform(SlowTrans)
        {
        }
      virtual ~DensityTransform()
        {
        }
      };

    //! Transform generalized model parameters to conductivity
    /*! Similarly to DensityTransform, this class transforms
     * generalized model parameters to conductivity taking an
     * intermediate step through slowness. This is motivated
     * by the parametrization used in the SINDRI project.
     */
    class ConductivityTransform: public jiba::GeneralModelTransform
      {
    private:
      const double a;
      const double b;
      const double c;
      boost::shared_ptr<GeneralModelTransform> SlownessTransform;
    public:
      //! Transform the normalized model parameters back to physical parameters
      virtual jiba::rvec GeneralizedToPhysical(const jiba::rvec &FullModel) const
        {
          jiba::rvec Slowness(SlownessTransform->GeneralizedToPhysical(
              FullModel));
          jiba::rvec Output(FullModel.size());
          for (size_t i = 0; i < FullModel.size(); ++i)
            {
              Output(i) = std::exp(-(a / (Slowness(i) * Slowness(i)) + b
                  / Slowness(i) + c));
            }
          return Output;
        }
      //! Transform the physical model parameters to generalized model parameters
      virtual jiba::rvec PhysicalToGeneralized(const jiba::rvec &FullModel) const
        {

          jiba::rvec Output(FullModel.size());
          for (size_t i = 0; i < FullModel.size(); ++i)
            {
              double vel = (-b + sqrt(b * b - 4.0 * a * (c + std::log(
                  FullModel(i))))) / (2 * a);
              Output(i) = 1.0 / vel;
            }
          return SlownessTransform->PhysicalToGeneralized(Output);
        }
      //! Transform the derivative with respect to the physical parameters to normalized parameters
      virtual jiba::rvec Derivative(const jiba::rvec &FullModel,
          const jiba::rvec &Derivative) const
        {
          jiba::rvec Slowness(SlownessTransform->GeneralizedToPhysical(
              FullModel));
          jiba::rvec SlowDeriv(SlownessTransform->Derivative(FullModel,
              Derivative));
          jiba::rvec Output(FullModel.size());
          for (size_t i = 0; i < FullModel.size(); ++i)
            {
              Output(i) = exp(-(a / (Slowness(i) * Slowness(i)) + b / Slowness(
                  i) + c)) * (2.0 * a * pow(Slowness(i), -3) + b / (Slowness(i)
                  * Slowness(i))) * SlowDeriv(i);
            }
          return Output;
        }
      //! The constructor needs a pointer to an object that gives slowness
      /*! To reduce the amount of code in the main program this class
       * takes a pointer to a transformation object that it uses to
       * transform the generalized model parameters to slowness before
       * then transforming to conductivity. In addition we can specify
       * the coefficients for the transformation.
       * @param SlowTrans A pointer to an object that gives slowness
       */
      ConductivityTransform(boost::shared_ptr<GeneralModelTransform> SlowTrans) :
        a(2.31e-7), b(-5.79e-4), c(0.124), SlownessTransform(VelTrans)
        {
        }
      virtual ~ConductivityTransform()
        {
        }
      };

    //!We can use this transformation class to chain a number of transformations together
    /*! This class takes a number of transformations and consecutively applies them
     * to the generalized parameters. For the backtransformation to physical
     * parameters the order is reversed and we use the chain rule to calculate
     * the derivatives.
     */
    class ChainedTransform: public jiba::GeneralModelTransform
      {
    private:
      std::vector<boost::shared_ptr<GeneralModelTransform> > Transforms;
    public:
      //! Transform the normalized model parameters back to physical parameters
      virtual jiba::rvec GeneralizedToPhysical(const jiba::rvec &FullModel) const
        {
          jiba::rvec Output(FullModel);

          for (size_t j = 0; j < Transforms.size(); ++j)
            Output = Transforms.at(j)->GeneralizedToPhysical(Output);
          return Output;
        }
      //! Transform the physical model parameters to generalized model parameters
      virtual jiba::rvec PhysicalToGeneralized(const jiba::rvec &FullModel) const
        {
          jiba::rvec Output(FullModel);
          for (int i = Transforms.size() - 1; i >= 0; --i)
            Output = Transforms.at(i)->PhysicalToGeneralized(Output);
          return Output;
        }
      //! Transform the derivative with respect to the physical parameters to normalized parameters
      virtual jiba::rvec Derivative(const jiba::rvec &FullModel,
          const jiba::rvec &Derivative) const
        {
          jiba::rvec Output(Derivative);
          jiba::rvec TransModel(FullModel);
          for (size_t j = 0; j < Transforms.size(); ++j)
            {
              Output = Transforms.at(j)->Derivative(TransModel, Output);
              TransModel = Transforms.at(j)->GeneralizedToPhysical(TransModel);
            }
          return Output;
        }
      //! Add a transform to the chain
      void AddTransform(boost::shared_ptr<GeneralModelTransform> Trans)
        {
          Transforms.push_back(Trans);
        }
      ChainedTransform()
        {

        }
      virtual ~ChainedTransform()
        {
        }
      };
  /* @} */
  }
#endif /* MODELTRANSFORMS_H_ */
