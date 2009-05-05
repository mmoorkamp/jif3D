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

class Rosenbrock: public jiba::ObjectiveFunction
  {
private:
  size_t neval;
  virtual void ImplDataDifference(const jiba::rvec &Model, jiba::rvec &Diff)
    {
      Diff.resize(2);
      double f1 = (Model(1) - Model(0) * Model(0));
      double f2 = 1. - Model(0);

      Diff(0) = 10 * f1;
      Diff(1) = f2;
      ++neval;
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
  size_t GetEval()
    {
      return neval;
    }
  Rosenbrock() :
    neval(0)
    {
    }
  virtual ~Rosenbrock()
    {
    }
  };

BOOST_AUTO_TEST_SUITE( OptStep_Test_Suite )

BOOST_AUTO_TEST_CASE (basic_nlcg_test)
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
          std::cout << std::endl;

        }
      const size_t neval = Objective->GetEval();
      BOOST_CHECK(neval < 100);
      BOOST_CHECK_CLOSE(Model(0),1.0,1.0);
      BOOST_CHECK_CLOSE(Model(1),1.0,1.0);
      std::cout << "Neval: " << neval << std::endl;
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
          std::cout << std::endl;

        }
      const size_t neval = Objective->GetEval();
      BOOST_CHECK(neval < 100);
      BOOST_CHECK_CLOSE(Model(0),1.0,0.01);
      BOOST_CHECK_CLOSE(Model(1),1.0,0.01);
      std::cout << "Neval: " << neval << std::endl;
    }

BOOST_AUTO_TEST_SUITE_END()
