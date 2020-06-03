//============================================================================
// Name        : test_inversion.h
// Author      : Max Moorkamp
// Version     :
// Copyright   : 2008, MM
//============================================================================

#define BOOST_TEST_MODULE Regularization test
#define BOOST_TEST_MAIN ...
#include "../Global/Jif3DTesting.h"
#include <boost/test/included/unit_test.hpp>
#include <stdlib.h>
#include <algorithm>
#include <numeric>
#include "MutualInformationConstraint.h"
#include <boost/test/tools/floating_point_comparison.hpp>

BOOST_AUTO_TEST_SUITE( Regularization_Test_Suite )

    void CheckGradient(jif3D::ObjectiveFunction &Objective, const jif3D::rvec &Model)
      {
        Objective.CalcMisfit(Model);
        jif3D::rvec Gradient = Objective.CalcGradient(Model);
        std::ofstream gradfile("migrad.out");
        for (size_t i = 0; i < Gradient.size(); ++i)
          {
            double delta = Model(i) * 0.001;
            jif3D::rvec Forward(Model);
            jif3D::rvec Backward(Model);
            Forward(i) += delta;
            Backward(i) -= delta;
            double FDGrad = (Objective.CalcMisfit(Forward)
                - Objective.CalcMisfit(Backward)) / (2 * delta);
            BOOST_CHECK_CLOSE(FDGrad, Gradient(i), 0.01);
            gradfile << i << " " << Gradient(i) << " " << FDGrad << " " << Model(i) << std::endl;

            /*if (std::abs(FDGrad) > 1e-10)
              {
                BOOST_CHECK_CLOSE(FDGrad, Gradient(i), 0.01);
              }
            else
              {
                BOOST_CHECK(std::abs(Gradient(i)) < 1e-10);
              }*/
          }
      }

    BOOST_AUTO_TEST_CASE (random_test)
      {
        srand48((unsigned int) time(nullptr));
        const size_t nx = 10;
        const size_t ny = 10;
        const size_t nz = 5;

        const size_t msize = nx * ny * nz * 2;
        jif3D::rvec StartModel(msize);
        std::generate(StartModel.begin(), StartModel.end(), drand48);

        jif3D::MutualInformationConstraint MI(0.0,1.0,0.0,1.0);

        double mi = MI.CalcMisfit(StartModel);
        std::cout << "Random MI: " << mi << std::endl;
        //jif3D::rvec Gradient = MI.CalcGradient(StartModel);
        //CheckGradient(MI, StartModel);

      }

    BOOST_AUTO_TEST_CASE (depend_test)
      {
        srand48((unsigned int) time(nullptr));
        const size_t nx = 10;
        const size_t ny = 10;
        const size_t nz = 5;

        const size_t msize = nx * ny * nz * 2;
        jif3D::rvec StartModel(msize);
        std::generate_n(StartModel.begin(), msize/2, drand48);
        for (size_t i = msize/2; i < msize; ++i)
          {
            StartModel(i) = StartModel(i-msize/2);
          }

        jif3D::MutualInformationConstraint MI(0.0,1.0,0.0,1.0);
        double mi = MI.CalcMisfit(StartModel);
        std::cout << "Depend MI: " << mi << std::endl;

        //CheckGradient(MI, StartModel);
      }

    BOOST_AUTO_TEST_CASE (depend_noise_test)
      {
        srand48((unsigned int) time(nullptr));
        const size_t nx = 7;
        const size_t ny = 4;
        const size_t nz = 3;

        const size_t msize = nx * ny * nz * 2;
        jif3D::rvec StartModel(msize);
        std::generate_n(StartModel.begin(), msize/2, drand48);
        for (size_t i = msize/2; i < msize; ++i)
          {
            StartModel(i) = StartModel(i-msize/2) + drand48();
          }

        jif3D::MutualInformationConstraint MI(0.0,1.0,0.0,2.0);
        double mi = MI.CalcMisfit(StartModel);
        std::cout << "Noise MI: " << mi << std::endl;

        //CheckGradient(MI, StartModel);
      }


    BOOST_AUTO_TEST_SUITE_END()
