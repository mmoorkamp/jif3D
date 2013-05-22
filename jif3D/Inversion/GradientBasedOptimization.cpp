//============================================================================
// Name        : GradientBasedOptimization.cpp
// Author      : Jun 29, 2009
// Version     : 
// Copyright   : 2009, mmoorkamp
//============================================================================


#include "GradientBasedOptimization.h"

namespace jif3D
  {

    GradientBasedOptimization::GradientBasedOptimization(boost::shared_ptr<
        jif3D::ObjectiveFunction> ObjFunction) :
      NonLinearOptimization(ObjFunction), HaveEvaluated_(false), RawGrad(),
          CovGrad(), SearchDir()
      {

      }

    GradientBasedOptimization::~GradientBasedOptimization()
      {

      }

    void GradientBasedOptimization::EvaluateModel(
        const jif3D::rvec &CurrentModel)
      {
        if (!HaveEvaluated_)
          {
            //these are very expensive so we only do it if we have to
            Misfit = GetObjective().CalcMisfit(CurrentModel);
            RawGrad = GetObjective().CalcGradient(CurrentModel);
          }
        CovGrad = ublas::element_prod(RawGrad, GetModelCovDiag());
        HaveEvaluated();
      }

  }
