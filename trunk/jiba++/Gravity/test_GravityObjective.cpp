//============================================================================
// Name        : test_GravityObjective.cpp
// Author      : Jun 6, 2009
// Version     :
// Copyright   : 2009, mmoorkamp
//============================================================================


#define BOOST_TEST_MODULE ThreeDGravityModel test
#define BOOST_TEST_MAIN ...
#include <boost/test/included/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/assign/std/vector.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/numeric/conversion/cast.hpp>

#include "test_common.h"
#include "GravityObjective.h"
#include "ThreeDGravityFactory.h"
#include "ThreeDGravityModel.h"
#include "FullSensitivityGravityCalculator.h"

BOOST_AUTO_TEST_SUITE( GravityObjective_Test_Suite )

void  CheckGradient(jiba::ObjectiveFunction &Objective, const jiba::rvec &Model)
    {
      Objective.CalcMisfit(Model);
      jiba::rvec Gradient = Objective.CalcGradient();

      BOOST_CHECK(Model.size() == Gradient.size());
      for (size_t i = 0; i < Gradient.size(); ++i)
        {
          double delta = Model(i) * 0.01;
          jiba::rvec Forward(Model);
          jiba::rvec Backward(Model);
          Forward(i) += delta;
          Backward(i) -= delta;
          if (delta> 0.0)
            {
              double FDGrad = (Objective.CalcMisfit(Forward) - Objective.CalcMisfit(Backward))/(2*delta);
              BOOST_CHECK_CLOSE(FDGrad,2*Gradient(i),0.001);
            }
        }
    }

  //check some general properties of the FTG tensor that hold for any measurement above the surface
  BOOST_AUTO_TEST_CASE (tensor_fdderiv_test)
    {
      jiba::ThreeDGravityModel GravityTest;
      const size_t nmeas = 3;
      const size_t ncells = 5;
      MakeRandomModel(GravityTest,ncells,nmeas,false);
      boost::shared_ptr<jiba::FullSensitivityGravityCalculator> TensorCalculator(jiba::CreateGravityCalculator<jiba::FullSensitivityGravityCalculator>::MakeTensor());
      jiba::rvec Observed(TensorCalculator->Calculate(GravityTest));
      Observed *= 1.1;

      boost::shared_ptr<jiba::GravityObjective> FTGObjective(
          new jiba::GravityObjective(true, false));
      FTGObjective->SetObservedData(Observed);
      FTGObjective->SetModelGeometry(GravityTest);
      //FTGObjective->SetDataCovar(jiba::ConstructError(FTGData, 0.02));
      //FTGObjective->SetPrecondDiag(PreCond);
      jiba::rvec InvModel(GravityTest.GetDensities().num_elements());
      std::copy(GravityTest.GetDensities().origin(),
          GravityTest.GetDensities().origin()
          + GravityTest.GetDensities().num_elements(), InvModel.begin());

      CheckGradient(*FTGObjective.get(),InvModel);
    }

  BOOST_AUTO_TEST_CASE (scalar_fdderiv_test)
    {
      jiba::ThreeDGravityModel GravityTest;
      const size_t nmeas = 3;
      const size_t ncells = 5;
      MakeRandomModel(GravityTest,ncells,nmeas,false);
      boost::shared_ptr<jiba::FullSensitivityGravityCalculator> TensorCalculator(jiba::CreateGravityCalculator<jiba::FullSensitivityGravityCalculator>::MakeScalar());
      jiba::rvec Observed(TensorCalculator->Calculate(GravityTest));
      Observed *= 1.1;

      boost::shared_ptr<jiba::GravityObjective> FTGObjective(
          new jiba::GravityObjective(false, false));
      FTGObjective->SetObservedData(Observed);
      FTGObjective->SetModelGeometry(GravityTest);
      //FTGObjective->SetDataCovar(jiba::ConstructError(FTGData, 0.02));
      //FTGObjective->SetPrecondDiag(PreCond);
      jiba::rvec InvModel(GravityTest.GetDensities().num_elements());
      std::copy(GravityTest.GetDensities().origin(),
          GravityTest.GetDensities().origin()
          + GravityTest.GetDensities().num_elements(), InvModel.begin());

      CheckGradient(*FTGObjective.get(),InvModel);
    }

  BOOST_AUTO_TEST_SUITE_END()
