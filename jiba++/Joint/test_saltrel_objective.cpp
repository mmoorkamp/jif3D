#define BOOST_TEST_MODULE SeismicModel test
#define BOOST_TEST_MAIN ...
#include <boost/test/included/unit_test.hpp>
#include <stdlib.h>
#include <time.h>
#include <vector>
#include "SaltRelConstraint.h"
#include "../Inversion/ModelTransforms.h"
#include "../Tomo/ThreeDSeismicModel.h"
#include "SaltRelTrans.h"
#include <boost/test/floating_point_comparison.hpp>
#include <boost/random/linear_congruential.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>

BOOST_AUTO_TEST_SUITE( SaltRel_Objective_Test_Suite )

jiba  ::rvec CheckGradient(jiba::ObjectiveFunction &Objective, const jiba::rvec &Model)
    {
      Objective.CalcMisfit(Model);
      jiba::rvec Gradient = Objective.CalcGradient(Model);
      jiba::rvec FDGrad(Model.size(),0.0);
      double Misfit = Objective.CalcMisfit(Model);
      for (size_t i = 0; i < Gradient.size(); ++i)
        {
          double delta = Model(i) * 0.0001;
          jiba::rvec Forward(Model);
          jiba::rvec Backward(Model);
          Forward(i) += delta;
          Backward(i) -= delta;
          FDGrad(i) = (Objective.CalcMisfit(Forward) - Objective.CalcMisfit(Backward))/(2.0 * delta);
          std::cout << i << std::endl;
          BOOST_CHECK_CLOSE(FDGrad(i),Gradient(i),0.1);
        }
      return FDGrad;
    }

  BOOST_AUTO_TEST_CASE (derivative_test)
    {
      const size_t ncells = 10;
      const size_t nparam = ncells * 3;
      jiba::rvec ModelVector(nparam);

  	jiba::ThreeDSeismicModel RelModel, ExclModel;
  	RelModel.SetSlownesses().resize(boost::extents[ncells][1][1]);
  	ExclModel.SetSlownesses().resize(boost::extents[ncells][1][1]);
  	std::fill_n(RelModel.SetSlownesses().origin(),ncells,1.0);
  	std::fill_n(ExclModel.SetSlownesses().origin(),ncells,0.0);

      boost::shared_ptr<jiba::GeneralModelTransform> DensTrans =
      boost::shared_ptr<jiba::GeneralModelTransform> (new jiba::DensityTransform(
              boost::shared_ptr<jiba::GeneralModelTransform>(new jiba::ModelCopyTransform()),RelModel));

      boost::shared_ptr<jiba::GeneralModelTransform> CondTrans =
      boost::shared_ptr<jiba::GeneralModelTransform> (new jiba::ConductivityTransform(
              boost::shared_ptr<jiba::GeneralModelTransform>(new jiba::ModelCopyTransform()),RelModel));

      const double minslow = 2e-4;
      const double maxslow = 5e-4;
      srand48(time(0));
      for (size_t i = 0; i < ncells; ++i)
        {
          ModelVector(i) = minslow + drand48() * (maxslow - minslow);
        }
      ublas::subrange(ModelVector,ncells,2*ncells) = DensTrans->GeneralizedToPhysical(ublas::subrange(ModelVector,0,ncells));
      ublas::subrange(ModelVector,2* ncells,3*ncells) = CondTrans->GeneralizedToPhysical(ublas::subrange(ModelVector,0,ncells));
      //check that constraint values are small when everything obeys the background relationship
      jiba::SaltRelConstraint SaltObjective(DensTrans,CondTrans);
      SaltObjective.SetExcludeCells(ExclModel);
      BOOST_CHECK_SMALL(SaltObjective.CalcMisfit(ModelVector),1e-12);
      //then also the gradient should be zero
      jiba::rvec ZeroGradient = SaltObjective.CalcGradient(ModelVector);
      for (size_t i = 0; i < nparam; ++i)
        {
          BOOST_CHECK_SMALL(ZeroGradient(i),1e-12);
        }
      srand(time(0));
      //now change one cell to salt values
      const size_t saltindex = rand() % ncells;
      ModelVector(saltindex) =1.0 / jiba::saltvel;
      ModelVector(saltindex + ncells) =jiba::saltdens;
      ModelVector(saltindex + 2 * ncells) =1.0 / jiba::saltres;
      //everything should still be effectively zero
      BOOST_CHECK_SMALL(SaltObjective.CalcMisfit(ModelVector),1e-12);
      //then also the gradient should be zero
      ZeroGradient = SaltObjective.CalcGradient(ModelVector);
      for (size_t i = 0; i < nparam; ++i)
        {
          BOOST_CHECK_SMALL(ZeroGradient(i),1e-12);
        }
      //now that we perturbed density and conductivity we should have a non-zero gradient
      subrange(ModelVector,ncells,3*ncells) *= 1.1;
      CheckGradient(SaltObjective,ModelVector);
    }

  BOOST_AUTO_TEST_SUITE_END()
