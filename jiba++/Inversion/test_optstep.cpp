//============================================================================
// Name        : test_optstep.cpp
// Author      : Apr 20, 2009
// Version     :
// Copyright   : 2009, mmoorkamp
//============================================================================


#define BOOST_TEST_MODULE OptStep test
#define BOOST_TEST_MAIN ...
#include <boost/test/included/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <stdlib.h>
#include <fstream>
#include "NonLinearConjugateGradient.h"
#include "LimitedMemoryQuasiNewton.h"
#include "ObjectiveFunction.h"
#include "JointObjective.h"
#include "../ModelTransforms/ModelCopyTransform.h"


class Rosenbrock: public jiba::ObjectiveFunction
  {
private:
  virtual void ImplDataDifference(const jiba::rvec &Model, jiba::rvec &Diff)
    {
      Diff.resize(2);
      double f1 = (Model(1) - Model(0) * Model(0));
      double f2 = 1. - Model(0);

      Diff(0) = 10 * f1;
      Diff(1) = f2;
    }
  virtual jiba::rvec ImplGradient(const jiba::rvec &Model,
      const jiba::rvec &Diff)
    {
      double f1 = (Model(1) - Model(0) * Model(0));
      double f2 = 1. - Model(0);

      jiba::rvec Gradient(2);
      Gradient(0) = -400. * f1 * Model(0) - 2. * f2;
      Gradient(1) = 200. * f1;
      return Gradient;
    }
public:
  Rosenbrock()
    {
    }
  virtual ~Rosenbrock()
    {
    }
  };

BOOST_AUTO_TEST_SUITE( OptStep_Test_Suite )

BOOST_AUTO_TEST_CASE  (basic_nlcg_test)
    {
      boost::shared_ptr<Rosenbrock> Objective(new Rosenbrock());
      jiba::NonLinearConjugateGradient NLCG(Objective);
      jiba::rvec Cov(2);
      //std::fill_n(Cov.begin(), 2, 1.0);
      Cov(0) = 10.0;
      Cov(1) = 0.1;
      NLCG.SetModelCovDiag(Cov);
      jiba::rvec Model(2);
      Model(0) = -1.2;
      Model(1) = 1.0;
      double Misfit = 1e10;
      while (Misfit > 1e-9)
        {
          NLCG.MakeStep(Model);
          Misfit = NLCG.GetMisfit();
        }
      const size_t neval = Objective->GetNEval();
      BOOST_CHECK(neval < 100);
      BOOST_CHECK_CLOSE(Model(0),1.0,0.01);
      BOOST_CHECK_CLOSE(Model(1),1.0,0.01);
      std::cout << "NLCG Neval: " << neval << std::endl;
    }

  BOOST_AUTO_TEST_CASE (basic_lbfgs_test)
    {
      boost::shared_ptr<Rosenbrock> Objective(new Rosenbrock());
      jiba::LimitedMemoryQuasiNewton LBFGS(Objective, 5);
      jiba::rvec Cov(2);
      //std::fill_n(Cov.begin(), 2, 1.0);
      Cov(0) = 10.0;
      Cov(1) = 0.1;
      LBFGS.SetModelCovDiag(Cov);

      jiba::rvec Model(2);
      Model(0) = -1.2;
      Model(1) = 1.0;
      double Misfit = 1e10;
      while (Misfit > 1e-9)
        {
          LBFGS.MakeStep(Model);
          Misfit = LBFGS.GetMisfit();
        }
      const size_t neval = Objective->GetNEval();
      BOOST_CHECK(neval < 100);
      BOOST_CHECK_CLOSE(Model(0),1.0,0.01);
      BOOST_CHECK_CLOSE(Model(1),1.0,0.01);
      std::cout << "LBFGS Neval: " << neval << std::endl;
    }

  BOOST_AUTO_TEST_CASE (basic_jointobjective_test)
    {
      boost::shared_ptr<jiba::JointObjective> Objective(new jiba::JointObjective());
      boost::shared_ptr<Rosenbrock> Rosen1(new Rosenbrock());
      boost::shared_ptr<Rosenbrock> Rosen2(new Rosenbrock());
      boost::shared_ptr<Rosenbrock> Rosen3(new Rosenbrock());

      boost::shared_ptr<jiba::GeneralModelTransform> Transform(new jiba::ModelCopyTransform());
      Objective->AddObjective(Rosen1,Transform,0.8);
      Objective->AddObjective(Rosen2,Transform,0.2);

      jiba::LimitedMemoryQuasiNewton JointLBFGS(Objective, 5);
      jiba::LimitedMemoryQuasiNewton RosenLBFGS(Rosen3, 5);
      jiba::rvec Cov(2);
      //std::fill_n(Cov.begin(), 2, 1.0);
      Cov(0) = 10.0;
      Cov(1) = 0.1;
      JointLBFGS.SetModelCovDiag(Cov);
      RosenLBFGS.SetModelCovDiag(Cov);

      jiba::rvec JointModel(2);
      JointModel(0) = -1.2;
      JointModel(1) = 1.0;
      jiba::rvec RosenModel(JointModel);
      double Misfit = 1e10;
      while (Misfit > 1e-9)
        {
          JointLBFGS.MakeStep(JointModel);
          RosenLBFGS.MakeStep(RosenModel);
          Misfit = RosenLBFGS.GetMisfit();
        }
      BOOST_CHECK_CLOSE(JointModel(0),RosenModel(0),0.01);
      BOOST_CHECK_CLOSE(JointModel(1),RosenModel(1),0.01);
    }
  BOOST_AUTO_TEST_SUITE_END()
