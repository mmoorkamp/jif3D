#define BOOST_TEST_MODULE SeismicModel test
#define BOOST_TEST_MAIN ...
#include <boost/test/included/unit_test.hpp>
#include <stdlib.h>
#include <time.h>
#include <vector>
#include "SaltRelConstraint.h"
#include "../Inversion/ModelTransforms.h"
#include "../Tomo/ThreeDSeismicModel.h"
#include "../Global/Jif3DPlatformHelper.h"
#include <boost/test/tools/floating_point_comparison.hpp>
#include <boost/random/linear_congruential.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>

BOOST_AUTO_TEST_SUITE( SaltRel_Objective_Test_Suite )

    jif3D::rvec CheckGradient(jif3D::ObjectiveFunction &Objective,
        const jif3D::rvec &Model)
      {
        Objective.CalcMisfit(Model);
        jif3D::rvec Gradient = Objective.CalcGradient(Model);
        jif3D::rvec FDGrad(Model.size(), 0.0);
        for (size_t i = 0; i < Gradient.size(); ++i)
          {
            double delta = Model(i) * 0.0001;
            jif3D::rvec Forward(Model);
            jif3D::rvec Backward(Model);
            Forward(i) += delta;
            Backward(i) -= delta;
            FDGrad(i) = (Objective.CalcMisfit(Forward) - Objective.CalcMisfit(Backward))
                / (2.0 * delta);
            BOOST_CHECK_CLOSE(FDGrad(i), Gradient(i), 0.1);
          }
        return FDGrad;
      }

    BOOST_AUTO_TEST_CASE (derivative_test)
      {
        const size_t ncells = 10;
        const size_t nparam = ncells * 3;
        jif3D::rvec ModelVector(nparam);

        jif3D::ThreeDSeismicModel RelModel, ExclModel;
        RelModel.SetSlownesses().resize(boost::extents[ncells][1][1]);
        ExclModel.SetSlownesses().resize(boost::extents[ncells][1][1]);
        std::fill_n(RelModel.SetSlownesses().origin(), ncells, 1.0);
        std::fill_n(ExclModel.SetSlownesses().origin(), ncells, 0.0);

        boost::shared_ptr<jif3D::GeneralModelTransform> DensTrans = boost::shared_ptr<
            jif3D::GeneralModelTransform>(
            new jif3D::DensityTransform(
                boost::shared_ptr<jif3D::GeneralModelTransform>(
                    new jif3D::ModelCopyTransform()), RelModel));

        boost::shared_ptr<jif3D::GeneralModelTransform> CondTrans = boost::shared_ptr<
            jif3D::GeneralModelTransform>(
            new jif3D::ConductivityTransform(
                boost::shared_ptr<jif3D::GeneralModelTransform>(
                    new jif3D::ModelCopyTransform()), RelModel));

        const double minslow = 2e-4;
        const double maxslow = 5e-4;
        jif3D::platform::srand48((unsigned int) time(nullptr));
        for (size_t i = 0; i < ncells; ++i)
          {
            ModelVector(i) = minslow + jif3D::platform::drand48() * (maxslow - minslow);
          }
        ublas::subrange(ModelVector, ncells, 2 * ncells) =
            DensTrans->GeneralizedToPhysical(ublas::subrange(ModelVector, 0, ncells));
        ublas::subrange(ModelVector, 2 * ncells, 3 * ncells) =
            CondTrans->GeneralizedToPhysical(ublas::subrange(ModelVector, 0, ncells));
        //check that constraint values are small when everything obeys the background relationship
        jif3D::SaltRelConstraint SaltObjective(DensTrans, CondTrans);
        SaltObjective.SetExcludeCells(ExclModel);
        BOOST_CHECK_SMALL(SaltObjective.CalcMisfit(ModelVector), 1e-12);
        //then also the gradient should be zero
        jif3D::rvec ZeroGradient = SaltObjective.CalcGradient(ModelVector);
        for (size_t i = 0; i < nparam; ++i)
          {
            BOOST_CHECK_SMALL(ZeroGradient(i), 1e-12);
          }
        srand((unsigned int) time(nullptr));
        //now change one cell to salt values
        const size_t saltindex = rand() % ncells;
        ModelVector(saltindex) = 1.0 / jif3D::saltvel;
        ModelVector(saltindex + ncells) = jif3D::saltdens;
        ModelVector(saltindex + 2 * ncells) = 1.0 / jif3D::saltres;
        //everything should still be effectively zero
        BOOST_CHECK_SMALL(SaltObjective.CalcMisfit(ModelVector), 1e-12);
        //then also the gradient should be zero
        ZeroGradient = SaltObjective.CalcGradient(ModelVector);
        for (size_t i = 0; i < nparam; ++i)
          {
            BOOST_CHECK_SMALL(ZeroGradient(i), 1e-12);
          }
        //now that we perturbed density and conductivity we should have a non-zero gradient
        subrange(ModelVector, ncells, 3 * ncells) *= 1.1;
        CheckGradient(SaltObjective, ModelVector);
      }

    BOOST_AUTO_TEST_SUITE_END()
