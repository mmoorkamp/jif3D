//============================================================================
// Name        : ThreeDModelObjective.h
// Author      : Nov 29, 2010
// Version     : 
// Copyright   : 2010, mmoorkamp
//============================================================================

#ifndef THREEDMODELOBJECTIVE_H_
#define THREEDMODELOBJECTIVE_H_

#include <boost/serialization/serialization.hpp>
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/export.hpp>
#include "../ModelBase/ModelRefiner.h"
#include "../Global/FatalException.h"
#include "../Inversion/ObjectiveFunction.h"

namespace jiba
  {
    //! This is a generic class for objective functions that use a Calculator object to calculate 3D synthetic data
    /*! The different types of synthetic data calculations we currently perform during the joint inversion
     * are all encapsulated in a Calculator object with identical interfaces. We therefore define this
     * template that calculates the misfit between observed and calculated data and provides other
     * common functionality with the type of the calculator object as a template parameter.
     * This template has 3 requirements on the calculator object.
     *   - A public typedef ModelType that specifies the model object type required for the forward calculation
     *   - A function  with signature void Calculate(const ModelType &); to calculate synthetic data for a given model
     *   - A function with signature jiba::rvec LQDerivative(const ModelType &,const jiba::rvec &); to calculate
     *      the derivative of a least squares objective function for a given model and misfit.
     * If all these requirements are met this class can use the calculator object to calculate misfit
     * and the derivative of the objective function.
     *
     */
    template<class ThreeDCalculatorType>
    class ThreeDModelObjective: public jiba::ObjectiveFunction
      {
    public:
      //! We create a shorthand for the template parameter that determines the type of the forward calculation object
      typedef ThreeDCalculatorType CalculatorType;
      //! The forward calculation class must contain a type definition that sets the type of the model object
      typedef typename CalculatorType::ModelType ModelType;
    protected:
      //! The object that calculates the synthetic data its type is set by the template parameter
      CalculatorType Calculator;
      //! The skeleton for the physical property model used in the inversion that contains geometry information etc.
      /*! The objective function simply copies the vector of values it receives for the current model
       * into this object, to obtain a well formed model for the synthetic forward and gradient calculation.
       * We therefore have to make sure that the total number of elements in the inversion grid
       * matches the size of the model vector. The type of this object is set by a nested type declaration
       * in the Calculator object.
       */
      ModelType CoarseModel;
      //! The geometry of the model used for forward calculations
      /*! If set by the user this object contains the geometry of a refined model
       *  that is used to calculate numerically correct synthetic data on the
       * coarse inversion grid. See also the Refiner object.
       */
      ModelType FineModel;
      //! Did we actually set a fine model for the forward calculation
      bool wantrefinement;
      //! We might use a different discretization for forward modeling and inversion
      /*! To ensure numerically correct results we might need a model discretization
       * that is much finer than the inversion grid we want to use. The Refiner object
       * takes CoarseModel which is the model with a geometry as intended in the
       * inversion and refines it using the geometry of FineModel. This model
       * is then used for the forward calculation.
       */
      jiba::ModelRefiner Refiner;
      //! The vector of observed data for all stations. The exact ordering of the data depends on the calculator object.
      jiba::rvec ObservedData;
      //! Calculate the difference between observed and synthetic data for a given model
      friend class boost::serialization::access;
      //! Provide serialization to be able to store objects and, more importantly for simpler MPI parallelization
      template<class Archive>
      void serialize(Archive & ar, const unsigned int version)
        {
          ar & boost::serialization::base_object<ObjectiveFunction>(*this);
          ar & Calculator;
          ar & CoarseModel;
          ar & FineModel;
          ar & wantrefinement;
          ar & Refiner;
          ar & ObservedData;
        }
    protected:
      ThreeDModelObjective()
        {

        }
    private:
      //The synthetic data from the last forward calculation
      jiba::rvec SynthData;
      virtual void
      ImplDataDifference(const jiba::rvec &Model, jiba::rvec &Diff);
      //! The implementation of the gradient calculation
      virtual jiba::rvec ImplGradient(const jiba::rvec &Model, const jiba::rvec &Diff);
      //! So far transformations have no effect
      virtual void SetDataTransformAction()
        {
          //TODO Implement transformation if necessary
        }
    public:
      //! Set the observed  data
      /*! We have to set the observed data in a format consistent with the calculator object.
       * As we provide an abstract interface here, the user has to make sure that the calculator
       * object yields a synthetic data vector with the same ordering as the observed data here, for
       * example by adding measurement points to the model object in the right order.
       * @param Data A vector with real values containing the observed data
       */
      void SetObservedData(const jiba::rvec &Data)
        {
          if (Data.empty())
            throw jiba::FatalException(
                "Cannot have empty observations in objective function.");
          ObservedData = Data;
        }
      //! Return a read only version of the observed data
      const jiba::rvec &GetObservedData() const
        {
          return ObservedData;
        }
      //! Return the synthetic data calculated for the last forward call
      /*! We do not store the synthetic data, but only the data difference which we need frequently.
       * As we only need the synthetic data occasionally, we generate it from the data difference
       * when needed and save memory.
       */
      jiba::rvec GetSyntheticData() const
        {
          return SynthData;
        }
      //! Set a skeleton for the  model that contains all information about cell sizes, site coordinates etc.
      /*! During the inversion we only copy a vector of values, so we have to store the
       * geometry and background information so we can form a complete model for the
       * forward calculation. The coarse model has to have the same number of cells as the model vector
       * that is passed to the gradient and forward calculations. It is used to establish the geometry
       * of the inversion model. For forward modeling this can be refined to the geometry of the FineConductivityModel
       * @param Model The skeleton object that contains all information about geometry, background and the location of the sites
       */
      void SetCoarseModelGeometry(const ModelType &Model)
        {
          if (Model.GetData().num_elements() == 0)
            throw jiba::FatalException("Cannot have empty model in objective function.");
          CoarseModel = Model;
        }
      //! Set the geometry of a finely discretized model to ensure numerical precision
      /*! The geometry we want to use for the inversion might not be fine enough
       * to ensure a proper numerical calculation. In these cases we can pass the geometry
       * of a refined model to the objective function. Note that the resulting model
       * will have the grid interfaces from both the coarse model and the fine model
       * as described for the ModelRefiner class. The resulting model must still have
       * equal grid sizes in each of the horizontal directions. This is not checked here.
       * @param Model An object containing the refined grid coordinates, conductivity values, site locations
       *        and back are ignored.
       */
      void SetFineModelGeometry(const ModelType &Model)
        {
          if (Model.GetData().num_elements() == 0)
            throw jiba::FatalException("Cannot have empty model in objective function.");
          FineModel = Model;
          //set the coordinates of the refiner object according to the fine model
          //we want to use for the forward calculation
          Refiner.SetXCoordinates(FineModel.GetXCoordinates());
          Refiner.SetYCoordinates(FineModel.GetYCoordinates());
          Refiner.SetZCoordinates(FineModel.GetZCoordinates());
          //set the flag that we call the refiner in the forward and gradient calculation
          wantrefinement = true;
        }

    public:
      //! The constructor needs a forward calculation object as a parameter
      /*! As this is a general objective function class for different types of forward calculations
       * we need to pass the object that calculates the synthetic data as
       * a parameter to the constructor. See the general description of
       * this class for requirements on the forward calculation object.
       * @param Calc The forward calculation object
       */
      explicit ThreeDModelObjective(const ThreeDCalculatorType &Calc);
      virtual ~ThreeDModelObjective();
      };

    template<class ThreeDCalculatorType>
    ThreeDModelObjective<ThreeDCalculatorType>::ThreeDModelObjective(
        const ThreeDCalculatorType &Calc) :
        Calculator(Calc), wantrefinement(false)
      {
      }

    template<class ThreeDCalculatorType>
    ThreeDModelObjective<ThreeDCalculatorType>::~ThreeDModelObjective()
      {
      }

    template<class ThreeDCalculatorType>
    void ThreeDModelObjective<ThreeDCalculatorType>::ImplDataDifference(
        const jiba::rvec &Model, jiba::rvec &Diff)
      {
        //make sure the sizes match
        assert(Model.size() == CoarseModel.GetData().num_elements());
        //Copy the model vector into the object with the geometry information
        std::copy(Model.begin(), Model.end(), CoarseModel.SetData().origin());

        //depending on whether we want to refine the inversion model
        //for the forward calculation or not we call the forward
        //calculation object with different models
        if (wantrefinement)
          {
            //project the values from the coarse model onto the fine model
            Refiner.RefineModel(CoarseModel, FineModel);
            //Calculate the synthetic data for the 3D model
            SynthData = Calculator.Calculate(FineModel);
          }
        else
          {
            SynthData = Calculator.Calculate(CoarseModel);
          }

        if (SynthData.size() != ObservedData.size())
          throw jiba::FatalException(
              " ThreeDModelObjective: Forward calculation does not give same amount of data !");
        Diff.resize(ObservedData.size());
        //calculate the difference between observed and synthetic
        std::transform(SynthData.begin(), SynthData.end(), ObservedData.begin(),
            Diff.begin(), std::minus<double>());
      }

    //The implementation of the gradient calculation
    template<class ThreeDCalculatorType>
    jiba::rvec ThreeDModelObjective<ThreeDCalculatorType>::ImplGradient(
        const jiba::rvec &Model, const jiba::rvec &Diff)
      {
        assert(Model.size() == CoarseModel.GetData().num_elements());
        //as we have stored the model vector from the misfit calculation
        //in the model object, we can check if the call is correct
        if (!std::equal(Model.begin(), Model.end(), CoarseModel.GetData().origin()))
          throw jiba::FatalException(
              "Gradient calculation needs identical model to forward !");
        //calculate the gradient
        if (wantrefinement)
          {
            Refiner.RefineModel(CoarseModel, FineModel);
            //calculate the gradient for the fine model
            jiba::rvec FineGradient(Calculator.LQDerivative(FineModel, Diff));
            //and return the projection of the fine gradient onto the coarse model
            return Refiner.CombineGradient(FineGradient, CoarseModel, FineModel);
          }
        //we only get here if we do not do any refinement
        //omitting the else saves us a compiler warning
        return Calculator.LQDerivative(CoarseModel, Diff);

      }

  }

#endif /* THREEDMODELOBJECTIVE_H_ */
