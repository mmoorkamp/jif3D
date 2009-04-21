//============================================================================
// Name        : GravityObjective.cpp
// Author      : Apr 16, 2009
// Version     :
// Copyright   : 2009, mmoorkamp
//============================================================================


#include "GravityObjective.h"
#include "MinMemGravityCalculator.h"
#include <cassert>
namespace jiba
  {

    GravityObjective::GravityObjective()
      {
        Calculator = boost::shared_ptr<jiba::ThreeDGravityCalculator>(
            jiba::CreateGravityCalculator<
                jiba::MinMemGravityCalculator>::MakeScalar());

      }

    GravityObjective::~GravityObjective()
      {

      }

    void GravityObjective::ImplDataDifference(const jiba::rvec &Model,
        jiba::rvec &Diff)
      {
        assert(DensityModel.GetMeasPosX().size() == ObservedData.size());
        std::copy(Model.begin(), Model.end(), DensityModel.SetDensities().origin());
        jiba::rvec SynthData(Calculator->Calculate(DensityModel));
        Diff.resize(ObservedData.size());
        std::transform(SynthData.begin(), SynthData.end(),
            ObservedData.begin(), Diff.begin(), std::minus<double>());
      }

    jiba::rvec GravityObjective::ImplGradient(const jiba::rvec &Model, const  jiba::rvec &Diff)
      {
        return Calculator->LQDerivative(DensityModel,Diff);
      }
  }
