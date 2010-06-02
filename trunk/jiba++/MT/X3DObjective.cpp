//============================================================================
// Name        : X3DObjective.cpp
// Author      : Jul 10, 2009
// Version     : 
// Copyright   : 2009, mmoorkamp
//============================================================================


#include "X3DObjective.h"
#include "../Global/FatalException.h"

namespace jiba
  {

    X3DObjective::X3DObjective() :
      wantrefinement(false)
      {

      }

    X3DObjective::~X3DObjective()
      {
      }

    void X3DObjective::ImplDataDifference(const jiba::rvec &Model,
        jiba::rvec &Diff)
      {
        assert(Model.size() == CoarseConductivityModel.GetConductivities().num_elements());
        //Copy the model vector into the object with the geometry information
        std::copy(Model.begin(), Model.end(),
            CoarseConductivityModel.SetConductivities().origin());
        //Calculate the impedances for the 3D model


        jiba::rvec SynthData;
        //project the values from the coarse model onto the fine model
        if (wantrefinement)
          {
            Refiner.RefineModel(CoarseConductivityModel, FineConductivityModel);
            //Calculate the travel times for the 3D model
            SynthData = Calculator.Calculate(FineConductivityModel);
          }
        else
          {
            SynthData = Calculator.Calculate(CoarseConductivityModel);
          }

        assert(SynthData.size() == ObservedData.size());
        Diff.resize(ObservedData.size());
        //calculate the difference between observed and synthetic
        std::transform(SynthData.begin(), SynthData.end(),
            ObservedData.begin(), Diff.begin(), std::minus<double>());
      }

    //The implementation of the gradient calculation
    jiba::rvec X3DObjective::ImplGradient(const jiba::rvec &Model,
        const jiba::rvec &Diff)
      {
        assert(Model.size() == CoarseConductivityModel.GetConductivities().num_elements());
        //as we have stored the model vector from the misfit calculation
        //in the model object, we can check if the call is correct
        if (!std::equal(Model.begin(), Model.end(),
            CoarseConductivityModel.GetConductivities().origin()))
          throw jiba::FatalException(
              "Gradient calculation needs identical model to forward !");
        //calculate the gradient

        if (wantrefinement)
          {
            Refiner.RefineModel(CoarseConductivityModel, FineConductivityModel);
            //calculate the gradient for the fine model
            jiba::rvec FineGradient(Calculator.LQDerivative(FineConductivityModel,
                Diff));
            //and return the projection of the fine gradient onto the coarse model
            return Refiner.CombineGradient(FineGradient, CoarseConductivityModel,
                FineConductivityModel);
          }
        //we only get here if we do not do any refinement
        //omitting the else saves us a compiler warning
        return Calculator.LQDerivative(CoarseConductivityModel, Diff);

      }
  }
