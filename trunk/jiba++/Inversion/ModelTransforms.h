//============================================================================
// Name        : ModelTransforms.h
// Author      : Apr 15, 2009
// Version     :
// Copyright   : 2009, mmoorkamp
//============================================================================

#ifndef MODELTRANSFORMS_H_
#define MODELTRANSFORMS_H_

#include <numeric>
#include <boost/serialization/serialization.hpp>
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/shared_ptr.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/math/special_functions/atanh.hpp>
#include "../Global/FatalException.h"
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
      //! The Reference model we devide the model parameters by
      const jiba::rvec Reference;
      friend class boost::serialization::access;
      //! Provide serialization to be able to store objects and, more importantly for simpler MPI parallelization
      template<class Archive>
      void serialize(Archive & ar, const unsigned int version)
        {
          ar & boost::serialization::base_object<GeneralModelTransform>(*this);
          ar & Reference;
        }
    public:
      //! We setup a clone function to have a virtual constructor and create polymorphic copies
      virtual NormalizeTransform* clone() const
        {
          return new NormalizeTransform(*this);
        }
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
          return GeneralizedToPhysical(Derivative);
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
      //! Each model parameter is divided by the reference values before taking the logarithm
      const jiba::rvec Reference;
      friend class boost::serialization::access;
      //! Provide serialization to be able to store objects and, more importantly for simpler MPI parallelization
      template<class Archive>
      void serialize(Archive & ar, const unsigned int version)
        {
          ar & boost::serialization::base_object<GeneralModelTransform>(*this);
          ar & Reference;
        }
    public:
      //! We setup a clone function to have a virtual constructor and create polymorphic copies
      virtual LogTransform* clone() const
        {
          return new LogTransform(*this);
        }
      //! Transform the normalized model parameters back to physical parameters
      virtual jiba::rvec GeneralizedToPhysical(const jiba::rvec &FullModel) const
        {
          assert(FullModel.size() == Reference.size());
          jiba::rvec Output(FullModel.size());
          for (size_t i = 0; i < FullModel.size(); ++i)
            Output(i) = std::exp(FullModel(i)) * Reference(i);
          return Output;
        }
      //! Transform the physical model parameters to generalized model parameters
      virtual jiba::rvec PhysicalToGeneralized(const jiba::rvec &FullModel) const
        {
          assert(FullModel.size() == Reference.size());
          jiba::rvec Output(FullModel.size());
          for (size_t i = 0; i < FullModel.size(); ++i)
            Output(i) = std::log(FullModel(i) / Reference(i));
          return Output;
        }
      //! Transform the derivative with respect to the physical parameters to normalized parameters
      virtual jiba::rvec Derivative(const jiba::rvec &FullModel,
          const jiba::rvec &Derivative) const
        {

          jiba::rvec Output(FullModel.size());
          for (size_t i = 0; i < FullModel.size(); ++i)
            {
              Output(i) = Reference(i) * std::exp(FullModel(i)) * Derivative(i);
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

    //! Change an unconstrained optimization problem to a constrained optimization problem through a tanh transformation
    /*! We often want to constrain the range of possible values in our optimization problem between an upper and
     * a lower limit. Instead of using a constrained optimization method, we use this transformation
     * \f$ m^{\star} =  \mbox{atanh} \left(2.0 * \frac{m - m_{min}} {m_{max} - m_{min}} \right) \f$ between the generalized
     * model parameters \f$ m \f$ and the generalized model parameters \f$ m^{\star} \f$. The generalized model
     * parameters can then vary over the whole range, while the physical model parameters will always been between
     * the two bounds. Note that with this transform the maximum and minimum values will be mapped to infinity
     * so in practice \f$ m_min < m < m_max \f$.
     */
    class TanhTransform: public jiba::GeneralModelTransform
      {
    private:
      //! The minimum value for the physical model parameters
      double min;
      //! The maximum value for the physical model parameters
      double max;
      friend class boost::serialization::access;
      //! Provide serialization to be able to store objects and, more importantly for simpler MPI parallelization
      template<class Archive>
      void serialize(Archive & ar, const unsigned int version)
        {
          ar & boost::serialization::base_object<GeneralModelTransform>(*this);
          ar & min;
          ar & max;
        }
    public:
      //! We setup a clone function to have a virtual constructor and create polymorphic copies
      virtual TanhTransform* clone() const
        {
          return new TanhTransform(*this);
        }
      //! Transform the normalized model parameters back to physical parameters
      virtual jiba::rvec GeneralizedToPhysical(const jiba::rvec &FullModel) const
        {
          jiba::rvec Output(FullModel.size());
          for (size_t i = 0; i < FullModel.size(); ++i)
            Output(i) = min + (1.0 + std::tanh(FullModel(i))) / 2.0 * (max - min);
          return Output;
        }
      //! Transform the physical model parameters to generalized model parameters
      virtual jiba::rvec PhysicalToGeneralized(const jiba::rvec &FullModel) const
        {
          jiba::rvec Output(FullModel.size());
          for (size_t i = 0; i < FullModel.size(); ++i)
            {
              if (FullModel(i) >= max || FullModel(i) <= min)
                {
                  std::cerr << i << " " << FullModel(i) << " " << max << " " << min
                      << std::endl;
                }
              const double argument = 2.0 * (FullModel(i) - min) / (max - min) - 1;
              Output(i) = boost::math::atanh(argument);
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
              Output(i) = (max - min) / (2.0 * std::pow(std::cosh(FullModel(i)), 2))
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
     * relationship between velocity and density
     * \f$ \rho = (1/s + b)/a \f$. This type of transformation
     * is motivated by the common velocity-density relationships
     * used in inversion problems and also observed for the Faeroe data.
     */
    class DensityTransform: public jiba::GeneralModelTransform
      {
    private:
      //! A pointer to a transformation that gives slowness
      boost::shared_ptr<GeneralModelTransform> SlownessTransform;
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
          ar & a;
          ar & b;
        }
      DensityTransform() :
          a(0), b(0)
        {
        }
    public:
      //! We setup a clone function to have a virtual constructor and create polymorphic copies
      virtual DensityTransform* clone() const
        {
          return new DensityTransform(*this);
        }
      //! Transform the normalized model parameters back to physical parameters, in this case from Slowness to Density
      virtual jiba::rvec GeneralizedToPhysical(const jiba::rvec &FullModel) const
        {
          jiba::rvec Slowness(SlownessTransform->GeneralizedToPhysical(FullModel));
          jiba::rvec Output(FullModel.size());
          for (size_t i = 0; i < FullModel.size(); ++i)
            {
              Output(i) = (1.0 / Slowness(i) + b) / a;
            }
          return Output;
        }
      //! Transform from Density to Slowness
      virtual jiba::rvec PhysicalToGeneralized(const jiba::rvec &FullModel) const
        {

          jiba::rvec Output(FullModel.size());
          for (size_t i = 0; i < FullModel.size(); ++i)
            {
              Output(i) = 1.0 / (a * FullModel(i) - b);
            }
          return SlownessTransform->PhysicalToGeneralized(Output);
        }
      //! Transform the derivative with respect to the Slowness to Density
      virtual jiba::rvec Derivative(const jiba::rvec &FullModel,
          const jiba::rvec &Derivative) const
        {
          jiba::rvec Slowness(SlownessTransform->GeneralizedToPhysical(FullModel));
          jiba::rvec SlowDeriv(SlownessTransform->Derivative(FullModel, Derivative));
          jiba::rvec Output(FullModel.size());
          for (size_t i = 0; i < FullModel.size(); ++i)
            {
              Output(i) = -1.0 / (Slowness(i) * Slowness(i)) / a * SlowDeriv(i);
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
       * @param aval The slope a of the functional relationship
       * @param bval The offset value (see equation above)
       */
      DensityTransform(boost::shared_ptr<GeneralModelTransform> SlowTrans, double aval =
          5000, double bval = 8500) :
          SlownessTransform(SlowTrans), a(aval), b(bval)
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
     * The equation for the transform is
     * \f$ \sigma = \exp( -a/s^2 - b/s -c), \f$
     * where s is slowness and the coefficient a,b,c can be specified
     * in the constructor.
     */
    class ConductivityTransform: public jiba::GeneralModelTransform
      {
    private:
      //! The coefficient for the quadratic term of the transform
      const double a;
      //! The coefficient for the linear term  of the transform
      const double b;
      //! The constant term of the transform
      const double c;
      boost::shared_ptr<GeneralModelTransform> SlownessTransform;
      friend class boost::serialization::access;
      //! Provide serialization to be able to store objects and, more importantly for simpler MPI parallelization
      template<class Archive>
      void serialize(Archive & ar, const unsigned int version)
        {
          ar & boost::serialization::base_object<GeneralModelTransform>(*this);
          ar & a;
          ar & b;
          ar & c;
          ar & SlownessTransform;
        }
    public:
      //! We setup a clone function to have a virtual constructor and create polymorphic copies
      virtual ConductivityTransform* clone() const
        {
          return new ConductivityTransform(*this);
        }
      //! Transform Generalized parameters in terms of slowness to conductivity using a functional relationship
      virtual jiba::rvec GeneralizedToPhysical(const jiba::rvec &FullModel) const
        {
          jiba::rvec Slowness(SlownessTransform->GeneralizedToPhysical(FullModel));
          jiba::rvec Conductivity(FullModel.size());
          for (size_t i = 0; i < FullModel.size(); ++i)
            {
              Conductivity(i) = std::exp(
                  -(a / (Slowness(i) * Slowness(i)) + b / Slowness(i) + c));
            }
          return Conductivity;
        }
      //! Transform Conductivity to Slowness and then Generalized Parameters
      virtual jiba::rvec PhysicalToGeneralized(const jiba::rvec &FullModel) const
        {
          size_t nvals = FullModel.size();
          jiba::rvec Slowness(nvals);
          for (size_t i = 0; i < nvals; ++i)
            {
              double vel = (-b + sqrt(b * b - 4.0 * a * (c + std::log(FullModel(i)))))
                  / (2.0 * a);
              Slowness(i) = 1.0 / vel;
            }
          return SlownessTransform->PhysicalToGeneralized(Slowness);
        }
      //! Transform the derivative with respect to the physical parameters to normalized parameters
      virtual jiba::rvec Derivative(const jiba::rvec &FullModel,
          const jiba::rvec &Derivative) const
        {
          jiba::rvec Slowness(SlownessTransform->GeneralizedToPhysical(FullModel));
          jiba::rvec SlowDeriv(SlownessTransform->Derivative(FullModel, Derivative));
          jiba::rvec Output(FullModel.size());
          for (size_t i = 0; i < FullModel.size(); ++i)
            {
              Output(i) =
                  std::exp(-(a / (Slowness(i) * Slowness(i)) + b / Slowness(i) + c))
                      * (2.0 * a * std::pow(Slowness(i), -3)
                          + b / (Slowness(i) * Slowness(i))) * SlowDeriv(i);
            }
          return Output;
        }
      //! The constructor needs a pointer to an object that gives slowness
      /*! To reduce the amount of code in the main program this class
       * takes a pointer to a transformation object that it uses to
       * transform the generalized model parameters to slowness before
       * then transforming to conductivity. In addition we can specify
       * the coefficients for the transformation.
       * We assume a functional relationship of the form
       * \f$ \rho = \exp( -a/v^2 + b/v +c) \f$.
       * @param SlowTrans A pointer to an object that gives slowness
       * @param aval The coefficient for the quadratic term
       * @param bval The coefficient for the linear term
       * @param cval The offset in the exponential
       */
      ConductivityTransform(boost::shared_ptr<GeneralModelTransform> SlowTrans,
          double aval = 2.31e-7, double bval = -5.79e-4, double cval = 0.124) :
          a(aval), b(bval), c(cval), SlownessTransform(SlowTrans)
        {
        }
      virtual ~ConductivityTransform()
        {
        }
      };

    //!We can use this transformation class to chain a number of transformations together
    /*! This class takes a number of transformations and consecutively applies them
     * to the generalized parameters. i.e. if we add the transforms f, g and h in that order,
     * we calculate the physical parameters p from the generalized parameters m as
     * \f$ p = h(g(f(m))) \f$ .
     * For the back transformation to physical
     * parameters the order is reversed and we use the chain rule to calculate
     * the derivatives.
     */
    class ChainedTransform: public jiba::GeneralModelTransform
      {
    private:
      //! We store pointers to each transform in the chain in this vector
      std::vector<boost::shared_ptr<GeneralModelTransform> > Transforms;
      friend class boost::serialization::access;
      //! Provide serialization to be able to store objects and, more importantly for simpler MPI parallelization
      template<class Archive>
      void serialize(Archive & ar, const unsigned int version)
        {
          ar & boost::serialization::base_object<GeneralModelTransform>(*this);
          ar & Transforms;
        }
    public:
      //! We setup a clone function to have a virtual constructor and create polymorphic copies
      virtual ChainedTransform* clone() const
        {
          return new ChainedTransform(*this);
        }
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
      //! Add a transform to the back of the chain
      void AppendTransform(boost::shared_ptr<GeneralModelTransform> Trans)
        {
          Transforms.push_back(Trans);
        }
      //! Add a transform to the front of the chain
      void PrependTransform(boost::shared_ptr<GeneralModelTransform> Trans)
        {
          Transforms.insert(Transforms.begin(), Trans);
        }
      ChainedTransform()
        {

        }
      virtual ~ChainedTransform()
        {
        }
      };


    //! For the cross-gradient we sometimes have to extract to sections from the model vector that are not continuous in memory
    /*! For cases were we want to extract several sections from the generalized model parameters
     * that are not continuous in memory, we cannot simply piece a transform
     * from two different SectionTransform objects. For example for the cross-gradient between seismic and MT
     * we currently store the model in the order slowness, density, conductivity. So we have to extract the
     * first n model parameters and the last n model parameters and form a vector of length 2n that
     * is fed to the cross-gradient function. This class generalizes this problem to an arbitrary number
     * of sections. In principle each section can have a different length n1, n2, ... and the resulting
     * transformed vector will have length n1+n2+...
     */
    class MultiSectionTransform: public jiba::GeneralModelTransform
      {
    private:
      size_t length;
      std::vector<size_t> startindices;
      std::vector<size_t> endindices;
      std::vector<boost::shared_ptr<GeneralModelTransform> > Transforms;
      friend class boost::serialization::access;
      //! Provide serialization to be able to store objects and, more importantly for simpler MPI parallelization
      template<class Archive>
      void serialize(Archive & ar, const unsigned int version)
        {
          ar & boost::serialization::base_object<GeneralModelTransform>(*this);
          ar & length;
          ar & startindices;
          ar & endindices;
          ar & Transforms;
        }
      MultiSectionTransform() :
          length(0)
        {

        }
    public:

      //! We setup a clone function to have a virtual constructor and create polymorphic copies
      virtual MultiSectionTransform* clone() const
        {
          return new MultiSectionTransform(*this);
        }
      virtual jiba::rvec GeneralizedToPhysical(const jiba::rvec &FullModel) const
        {
          using boost::numeric::ublas::subrange;
    	  //we calculate the total length of the transformed vector from
    	  //the length of each section
          const size_t outlength = std::accumulate(endindices.begin(), endindices.end(),
              0) - std::accumulate(startindices.begin(), startindices.end(), 0);
          jiba::rvec Result(outlength);

          const size_t nsections = startindices.size();
          //the result will be continuous in memory and always
          //start at 0, each section will just follow after the other
          size_t resultstartindex = 0;
          for (size_t i = 0; i < nsections; ++i)
            {
              const size_t resultendindex = resultstartindex + endindices[i]
                  - startindices[i];
              subrange(Result, resultstartindex, resultendindex) =
                  Transforms[i]->GeneralizedToPhysical(
                      subrange(FullModel, startindices[i], endindices[i]));
              resultstartindex = resultendindex;
            }
          return Result;
        }
      virtual jiba::rvec PhysicalToGeneralized(const jiba::rvec &FullModel) const
        {
    	  // !! THIS NEEDS TO BE CHECKED
          return FullModel;
        }
      virtual jiba::rvec Derivative(const jiba::rvec &FullModel,
          const jiba::rvec &Derivative) const
        {
          using boost::numeric::ublas::subrange;
          jiba::rvec Result(length);
          Result.clear();
          const size_t nsections = startindices.size();
          size_t currstart = 0;
          for (size_t i = 0; i < nsections; ++i)
            {
              subrange(Result, startindices[i], endindices[i]) =
                  Transforms[i]->Derivative(
                      subrange(FullModel, startindices[i], endindices[i]),
                      subrange(Derivative, currstart,
                          currstart + endindices[i] - startindices[i]));
              currstart += endindices[i] - startindices[i];
            }
          return Result;
        }
      void AddSection(size_t startindex, size_t endindex,
          boost::shared_ptr<GeneralModelTransform> Trans)
        {
          startindices.push_back(startindex);
          endindices.push_back(endindex);
          Transforms.push_back(Trans);
        }
      MultiSectionTransform(size_t l) :
          length(l)
        {

        }
      MultiSectionTransform(size_t l, size_t startindex, size_t endindex,
              boost::shared_ptr<GeneralModelTransform> Trans) :
               length(l)
      {
    	  AddSection(startindex,endindex,Trans);
      }
      virtual ~MultiSectionTransform()
        {

        }
      };


  /* @} */
  }
#endif /* MODELTRANSFORMS_H_ */
