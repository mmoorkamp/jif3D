//============================================================================
// Name        : test_inversion.h
// Author      : Max Moorkamp
// Version     :
// Copyright   : 2008, MM
//============================================================================

#define BOOST_TEST_MODULE Inversion test
#define BOOST_TEST_MAIN ...
#include <boost/test/included/unit_test.hpp>
#include <stdlib.h>
#include <algorithm>
#include <numeric>
#include "../Gravity/test_common.h"
#include "MinDiffRegularization.h"
#include "GradientRegularization.h"
#include <boost/test/floating_point_comparison.hpp>

BOOST_AUTO_TEST_SUITE( Regularization_Test_Suite )

void  CheckGradient(jiba::ObjectiveFunction &Objective, const jiba::rvec &Model)
    {
      Objective.CalcMisfit(Model);
      jiba::rvec Gradient = Objective.CalcGradient();
      for (size_t i = 0; i < Gradient.size(); ++i)
        {
          double delta = Model(i) * 0.01;
          jiba::rvec Forward(Model);
          jiba::rvec Backward(Model);
          Forward(i) += delta;
          Backward(i) -= delta;
          double FDGrad = (Objective.CalcMisfit(Forward) - Objective.CalcMisfit(Backward))/(2*delta);
          BOOST_CHECK_CLOSE(FDGrad,Gradient(i),0.001);
        }
    }

  BOOST_AUTO_TEST_CASE (mindiff_test)
    {
      srand(time(NULL));
      const size_t msize = 18;
      jiba::rvec StartModel(msize), PertModel(msize);
      jiba::rvec PreCond(msize);
      std::fill(PreCond.begin(),PreCond.end(),1.0);
      std::generate(StartModel.begin(),StartModel.end(),rand);
      std::generate(PertModel.begin(),PertModel.end(),rand);

      jiba::MinDiffRegularization Regularization;
      Regularization.SetPrecondDiag(PreCond);
      Regularization.SetReferenceModel(StartModel);
      jiba::rvec Diff = StartModel - PertModel;
      double Misfit = Regularization.CalcMisfit(PertModel);
      BOOST_CHECK_CLOSE(Misfit,ublas::inner_prod(Diff,Diff),0.001);
      CheckGradient(Regularization,PertModel);
    }

  BOOST_AUTO_TEST_CASE (gradreg_test)
    {
      srand(time(NULL));
      jiba::ThreeDGravityModel GravModel;
      GravModel.SetDensities().resize(boost::extents[5][4][3]);

      const size_t msize = GravModel.GetDensities().num_elements();
      jiba::rvec StartModel(msize), PertModel(msize);
      jiba::rvec PreCond(msize);
      std::fill(PreCond.begin(),PreCond.end(),1.0);
      std::generate(StartModel.begin(),StartModel.end(),rand);
      std::generate(PertModel.begin(),PertModel.end(),rand);

      jiba::GradientRegularization Regularization(GravModel);
      Regularization.SetPrecondDiag(PreCond);
      Regularization.SetReferenceModel(StartModel);

      jiba::rvec Diff = StartModel - PertModel;
      double Misfit = Regularization.CalcMisfit(PertModel);

      double zero = Regularization.CalcMisfit(StartModel+PreCond);
      BOOST_CHECK_CLOSE(zero,0.0,0.0001);
      CheckGradient(Regularization,PertModel);
    }

  BOOST_AUTO_TEST_SUITE_END()
