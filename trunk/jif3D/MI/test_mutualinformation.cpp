//============================================================================
// Name        : test_inversion.h
// Author      : Max Moorkamp
// Version     :
// Copyright   : 2008, MM
//============================================================================

#define BOOST_TEST_MODULE MutualInformation test
#define BOOST_TEST_MAIN ...
#include "../Global/Jif3DTesting.h"
#include <boost/test/included/unit_test.hpp>
#include <boost/test/tools/floating_point_comparison.hpp>

#include <stdlib.h>
#include <algorithm>
#include <numeric>
#include "MutualInformationConstraint.h"

BOOST_AUTO_TEST_SUITE (MutualInformation_Test_Suite)

void CheckGradient(jif3D::ObjectiveFunction &Objective, const jif3D::rvec &Model)
  {
    Objective.CalcMisfit(Model);
    jif3D::rvec Gradient = Objective.CalcGradient(Model);
    std::ofstream gradfile("migrad.out");
    for (size_t i = 0; i < Gradient.size(); ++i)
      {
        double delta = 0.0001;
        jif3D::rvec Forward(Model);
        jif3D::rvec Backward(Model);
        Forward(i) += delta;
        Backward(i) -= delta;
        double FDGrad = (Objective.CalcMisfit(Forward) - Objective.CalcMisfit(Backward))
            / (2 * delta);
        BOOST_CHECK_CLOSE(FDGrad, Gradient(i), 1.0);
        gradfile << i << " " << Gradient(i) << " " << FDGrad << " " << Model(i)
            << std::endl;

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

BOOST_AUTO_TEST_CASE (se_diff_test)
  {
    srand48((unsigned int) time(nullptr));
    const size_t nx = 20;
    jif3D::rvec x(nx);
    std::generate(x.begin(),x.end(),drand48);

    double H = jif3D::shan_entropy(x);
    jif3D::rvec dx = jif3D::diff_shan_entropy(x);
    const double delta = 0.001;
    for (size_t i = 0; i < nx; ++i)
      {
        jif3D::rvec xtest(x);
        xtest(i) += delta;
        double pdx = jif3D::shan_entropy(xtest);
        xtest(i) -= 2*delta;
        double ndx = jif3D::shan_entropy(xtest);
        double fd = (pdx - ndx)/(2.0 * delta);
        BOOST_CHECK_CLOSE(fd,dx(i),0.1);
      }

  }

BOOST_AUTO_TEST_CASE (random_test)
  {
    srand48((unsigned int) time(nullptr));
    const size_t nx = 5;
    const size_t ny = 4;
    const size_t nz = 3;

    const size_t msize = nx * ny * nz * 2;
    jif3D::rvec StartModel(msize);
    std::generate(StartModel.begin(), StartModel.end(), drand48);

    jif3D::MutualInformationConstraint MI(0.0, 1.0, 0.0, 1.0, 50);

    double mi = MI.CalcMisfit(StartModel);
    std::cout << "Random MI: " << mi << std::endl;
    //jif3D::rvec Gradient = MI.CalcGradient(StartModel);
    CheckGradient(MI, StartModel);

  }

BOOST_AUTO_TEST_CASE (depend_test)
  {
    srand48((unsigned int) time(nullptr));
    const size_t nx = 5;
    const size_t ny = 4;
    const size_t nz = 3;

    const size_t msize = nx * ny * nz * 2;
    jif3D::rvec StartModel(msize);
    std::generate_n(StartModel.begin(), msize/2, drand48);
    for (size_t i = msize/2; i < msize; ++i)
      {
        StartModel(i) = StartModel(i-msize/2);
      }

    jif3D::MutualInformationConstraint MI(0.0, 1.0, 0.0, 1.0, 50);
    double mi = MI.CalcMisfit(StartModel);
    std::cout << "Depend MI: " << mi << std::endl;

    CheckGradient(MI, StartModel);
  }

BOOST_AUTO_TEST_CASE (depend_noise_test)
  {
    srand48((unsigned int) time(nullptr));
    const size_t nx = 5;
    const size_t ny = 4;
    const size_t nz = 3;

    const size_t msize = nx * ny * nz * 2;
    jif3D::rvec StartModel(msize);
    std::generate_n(StartModel.begin(), msize/2, drand48);
    for (size_t i = msize/2; i < msize; ++i)
      {
        StartModel(i) = StartModel(i-msize/2) + 0.1 * drand48();
      }

    jif3D::MutualInformationConstraint MI(0.0, 1.0, 0.0, 2.0, 50);
    double mi = MI.CalcMisfit(StartModel);
    std::cout << "Noise MI: " << mi << std::endl;

    CheckGradient(MI, StartModel);
  }

BOOST_AUTO_TEST_SUITE_END()
