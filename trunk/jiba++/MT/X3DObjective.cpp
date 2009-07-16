//============================================================================
// Name        : X3DObjective.cpp
// Author      : Jul 10, 2009
// Version     : 
// Copyright   : 2009, mmoorkamp
//============================================================================


#include "X3DObjective.h"
#include "X3DMTCalculator.h"

namespace jiba
  {

    X3DObjective::X3DObjective()
      {

      }

    X3DObjective::~X3DObjective()
      {
      }

    void X3DObjective::ImplDataDifference(const jiba::rvec &Model,
        jiba::rvec &Diff)
      {
        assert(Model.size() == ConductivityModel.GetConductivities().num_elements());
        //Copy the model vector into the object with the geometry information
        std::copy(Model.begin(), Model.end(),
            ConductivityModel.SetConductivities().origin());
        //Calculate the travel times for the 3D model
        jiba::X3DMTCalculator Calculator;
        jiba::rvec SynthData(Calculator.Calculate(ConductivityModel));
        Diff.resize(ObservedData.size());
        assert(SynthData.size() == ObservedData.size());
        //calculate the difference between observed and synthetic
        std::transform(SynthData.begin(), SynthData.end(),
            ObservedData.begin(), Diff.begin(), std::minus<double>());
      }

    //The implementation of the gradient calculation
    jiba::rvec X3DObjective::ImplGradient(const jiba::rvec &Model,
        const jiba::rvec &Diff)
      {
        assert(Model.size() == ConductivityModel.GetConductivities().num_elements());
        //Copy the model vector into the object with the geometry information
        std::copy(Model.begin(), Model.end(),
            ConductivityModel.SetConductivities().origin());
        //calculate the gradient
        jiba::X3DMTCalculator Calculator;
        return Calculator.LQDerivative(ConductivityModel, Diff);
      }
  }
