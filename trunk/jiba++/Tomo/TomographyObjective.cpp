//============================================================================
// Name        : TomographyObjective.cpp
// Author      : May 5, 2009
// Version     :
// Copyright   : 2009, mmoorkamp
//============================================================================


#include "TomographyObjective.h"

namespace jiba
  {

    TomographyObjective::TomographyObjective() :
      SlownessModel(), ObservedData(), Calculator()
      {
      }

    TomographyObjective::~TomographyObjective()
      {
      }

    void TomographyObjective::ImplDataDifference(const jiba::rvec &Model,
        jiba::rvec &Diff)
      {
        std::copy(Model.begin(), Model.end(),
            SlownessModel.SetSlownesses().origin());
        jiba::rvec SynthData(Calculator.Calculate(SlownessModel));
        Diff.resize(ObservedData.size());
        std::transform(SynthData.begin(), SynthData.end(),
            ObservedData.begin(), Diff.begin(), std::minus<double>());
      }

    jiba::rvec TomographyObjective::ImplGradient(const jiba::rvec &Model,
        const jiba::rvec &Diff)
      {
        std::copy(Model.begin(), Model.end(),
            SlownessModel.SetSlownesses().origin());
        return Calculator.LQDerivative(SlownessModel, Diff);
      }
  }
