//============================================================================
// Name        : ThreeDModelObjective.h
// Author      : Nov 29, 2010
// Version     :
// Copyright   : 2010, mmoorkamp
//============================================================================

#ifndef THREEDMODELOBJECTIVE_H_
#define THREEDMODELOBJECTIVE_H_

#include <fstream>

#include "../Global/Serialization.h"
#include "../Global/Jif3DGlobal.h"
#include "../Global/FatalException.h"
#include "../Inversion/ObjectiveFunction.h"
#include "OneDMTCalculator.h"

namespace jif3D
  {
    //! Calculate the data misfit for MT observations assuming a 1D layered half-space as Earth structure.
    class J3DEXPORT OneDMTObjective: public jif3D::ObjectiveFunction
      {

    private:
      jif3D::OneDMTCalculator Calculator;
      jif3D::X3DModel MTModel;
      jif3D::rvec ObservedData;
      //! Calculate the difference between observed and synthetic data for a given model
      friend class access;
      //! Provide serialization to be able to store objects and, more importantly for simpler MPI parallelization
      template<class Archive>
      void serialize(Archive & ar, const unsigned int version)
        {
          ar & base_object<ObjectiveFunction>(*this);
          ar & Calculator;
          ar & ObservedData;
        }
    private:
      virtual void ImplDataDifference(const jif3D::rvec &Model, jif3D::rvec &Diff) override
        {
          //make sure the sizes match
          //std::cout << "Model: " << Model << std::endl;
          assert(Model.size() == MTModel.GetBackgroundConductivities().size());
          //Copy the model vector into the object with the geometry information
          MTModel.SetBackgroundConductivities(
              std::vector<double>(Model.begin(), Model.end()));
          jif3D::rvec SynthData;
          SynthData = Calculator.Calculate(MTModel);
          if (SynthData.size() != ObservedData.size())
            throw jif3D::FatalException(
                " ThreeDModelObjective: Forward calculation does not give same amount of data !",
                __FILE__, __LINE__);
          Diff.resize(ObservedData.size());
          //calculate the difference between observed and synthetic
          std::transform(SynthData.begin(), SynthData.end(), ObservedData.begin(),
              Diff.begin(), std::minus<double>());
        }
      //! The implementation of the gradient calculation
      virtual jif3D::rvec ImplGradient(const jif3D::rvec &Model, const jif3D::rvec &Diff) override
        {
          double misfit = CalcMisfit(Model);
          const size_t nlayers = Model.size();
          jif3D::rvec result(nlayers);
          for (size_t index = 0; index < nlayers; ++index)
            {
              double delta = Model(index) * 0.01;
              jif3D::rvec Forward(Model);
              jif3D::rvec Backward(Model);
              Forward(index) += delta;
              Backward(index) -= delta;
              double ForFDGrad = (CalcMisfit(Forward) - misfit) / (delta);
              double BackFDGrad = (misfit - CalcMisfit(Backward)) / delta;
              result(index) = (ForFDGrad + BackFDGrad) / 2.0;

            }
          return result;
          assert(Model.size() == MTModel.GetBackgroundConductivities().size());
          //as we have stored the model vector from the misfit calculation
          //in the model object, we can check if the call is correct

          if (!std::equal(Model.begin(), Model.end(),
              MTModel.GetBackgroundConductivities().begin()))
            throw jif3D::FatalException(
                "Gradient calculation needs identical model to forward !", __FILE__,
                __LINE__);
          //calculate the gradient
          jif3D::rvec Deriv = Calculator.LQDerivative(MTModel, Diff);
          //std::cout << Calculator.GetGamma() << std::endl;
          //std::cout << "Deriv: " << Deriv << std::endl;

          return Deriv;
        }

    public:
      //! The clone function provides a virtual constructor
      virtual OneDMTObjective *clone() const override
        {
          return new OneDMTObjective(*this);
        }
      //! Set the observed  data
      /*! We have to set the observed data in a format consistent with the calculator object.
       * As we provide an abstract interface here, the user has to make sure that the calculator
       * object yields a synthetic data vector with the same ordering as the observed data here, for
       * example by adding measurement points to the model object in the right order.
       * @param Data A vector with real values containing the observed data
       */
      void SetObservedData(const jif3D::rvec &Data)
        {
          if (Data.empty())
            throw jif3D::FatalException(
                "Cannot have empty observations in objective function.", __FILE__,
                __LINE__);
          ObservedData = Data;
        }
      //! Return a read only version of the observed data
      const jif3D::rvec &GetObservedData() const
        {
          return ObservedData;
        }
      //! Return the synthetic data calculated for the last forward call
      /*! We do not store the synthetic data, but only the data difference which we need frequently.
       * As we only need the synthetic data occasionally, we generate it from the data difference
       * when needed and save memory.
       */
      jif3D::rvec GetSyntheticData() const
        {
          jif3D::rvec Synthetic(GetDataDifference());
          Synthetic = ublas::element_prod(Synthetic, GetDataError());
          Synthetic += ObservedData;
          return Synthetic;
        }
      //! Set a skeleton for the  model that contains all information about cell sizes, site coordinates etc.
      /*! During the inversion we only copy a vector of values, so we have to store the
       * geometry and background information so we can form a complete model for the
       * forward calculation. The coarse model has to have the same number of cells as the model vector
       * that is passed to the gradient and forward calculations. It is used to establish the geometry
       * of the inversion model. For forward modeling this can be refined to the geometry of the FineConductivityModel
       * @param Model The skeleton object that contains all information about geometry, background and the location of the sites
       */
      void SetModelGeometry(const jif3D::X3DModel &TheModel)
        {
          if (TheModel.GetBackgroundConductivities().size() == 0)
            throw jif3D::FatalException("Cannot have empty model in objective function.",
                __FILE__, __LINE__);
          MTModel = TheModel;
        }
    public:
      //! The constructor needs a forward calculation object as a parameter
      /*! As this is a general objective function class for different types of forward calculations
       * we need to pass the object that calculates the synthetic data as
       * a parameter to the constructor. See the general description of
       * this class for requirements on the forward calculation object.
       * @param Calc The forward calculation object
       */
      OneDMTObjective()
        {
        }
      virtual ~OneDMTObjective()
        {
        }
      };

  }

#endif /* THREEDMODELOBJECTIVE_H_ */
