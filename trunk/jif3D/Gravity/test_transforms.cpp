//============================================================================
// Name        : test_transforms.cpp
// Author      : Apr 22, 2009
// Version     :
// Copyright   : 2009, mmoorkamp
//============================================================================

#define BOOST_TEST_MODULE GravityTransforms test
#define BOOST_TEST_MAIN ...
#include "../Global/Jif3DTesting.h"
#include <boost/test/included/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <stdlib.h>
#include <time.h>
#include "test_common.h"
#include "GravityTransforms.h"
#include "../GravMag/FullSensitivityGravMagCalculator.h"
#include "ThreeDGravityFactory.h"
#include "../GravMag/MinMemGravMagCalculator.h"
#include "../Global/NumUtil.h"

BOOST_AUTO_TEST_SUITE( Gravity_Transforms_Suite )

    BOOST_AUTO_TEST_CASE(ftginvar_test)
      {
        jif3D::FTGInvariant InvarTrans;
        jif3D::rvec TestInput(9);
        TestInput(0) = 1.0;
        TestInput(1) = 2.0;
        TestInput(2) = 3.0;
        TestInput(3) = 2.0;
        TestInput(4) = 4.0;
        TestInput(5) = 5.0;
        TestInput(6) = 3.0;
        TestInput(7) = 5.0;
        TestInput(8) = 6.0;
        jif3D::rvec Invar(InvarTrans.Transform(TestInput));
        jif3D::rmat Deriv(InvarTrans.Derivative(TestInput));
        BOOST_CHECK(Invar.size() == 1);
        BOOST_CHECK(Deriv.size1() == 1);
        BOOST_CHECK(Deriv.size2() == 9);
        BOOST_CHECK_CLOSE(Invar(0), -4.0, std::numeric_limits<float>::epsilon());
        BOOST_CHECK_CLOSE(Deriv(0, 0), 10.0, std::numeric_limits<float>::epsilon());
        BOOST_CHECK_CLOSE(Deriv(0, 1), -2.0, std::numeric_limits<float>::epsilon());
        BOOST_CHECK_CLOSE(Deriv(0, 2), -3.0, std::numeric_limits<float>::epsilon());
        BOOST_CHECK_CLOSE(Deriv(0, 3), -2.0, std::numeric_limits<float>::epsilon());
        BOOST_CHECK_CLOSE(Deriv(0, 4), 7.0, std::numeric_limits<float>::epsilon());
        BOOST_CHECK_CLOSE(Deriv(0, 5), -5.0, std::numeric_limits<float>::epsilon());
        BOOST_CHECK_CLOSE(Deriv(0, 6), -3.0, std::numeric_limits<float>::epsilon());
        BOOST_CHECK_CLOSE(Deriv(0, 7), -5.0, std::numeric_limits<float>::epsilon());
        BOOST_CHECK_CLOSE(Deriv(0, 8), 5.0, std::numeric_limits<float>::epsilon());

      }

    BOOST_AUTO_TEST_CASE(ftginvar_sens_test)
      {
        typedef typename jif3D::FullSensitivityGravMagCalculator<jif3D::ThreeDGravityModel> CalculatorType;
        boost::shared_ptr<CalculatorType> InvCalculator =
            boost::shared_ptr<CalculatorType>(
                jif3D::CreateGravityCalculator<CalculatorType>::MakeTensor());
        boost::shared_ptr<jif3D::VectorTransform> Transform = boost::make_shared<jif3D::FTGInvariant>();
        InvCalculator->SetDataTransform(Transform);

        jif3D::ThreeDGravityModel GravityTest;

        const size_t nmeas = 10;
        const size_t ncells = 10;
        MakeRandomModel(GravityTest, ncells, nmeas, true);

        jif3D::rvec InvDirect(InvCalculator->Calculate(GravityTest));
        jif3D::rvec InvCached(InvCalculator->Calculate(GravityTest));
        BOOST_CHECK(InvDirect.size() == InvCached.size());
        for (size_t i = 0; i < InvDirect.size(); ++i)
          {
            BOOST_CHECK_CLOSE(InvDirect(i), InvCached(i),
                std::numeric_limits<float>::epsilon());
          }

      }

    BOOST_AUTO_TEST_CASE(ftginvar_deriv_test)
      {
        typedef typename jif3D::FullSensitivityGravMagCalculator<jif3D::ThreeDGravityModel> CalculatorType;
        boost::shared_ptr<CalculatorType> CacheCalculator = boost::shared_ptr<
            CalculatorType>(jif3D::CreateGravityCalculator<CalculatorType>::MakeTensor());

        typedef typename jif3D::FullSensitivityGravMagCalculator<jif3D::ThreeDGravityModel> CompType;
        boost::shared_ptr<CompType> CompCalculator = boost::shared_ptr<CompType>(
            jif3D::CreateGravityCalculator<CompType>::MakeTensor());

        boost::shared_ptr<jif3D::VectorTransform> Transform = boost::make_shared<
            jif3D::FTGInvariant>();
        CacheCalculator->SetDataTransform(Transform);
        CompCalculator->SetDataTransform(Transform);

        jif3D::ThreeDGravityModel GravityTest;

        const size_t nmeas = 10;
        const size_t ncells = 10;
        MakeRandomModel(GravityTest, ncells, nmeas, true);

        jif3D::rvec Misfit(nmeas);
        std::fill(Misfit.begin(), Misfit.end(), 1.0);
        jif3D::rvec InvDirect(CompCalculator->LQDerivative(GravityTest, Misfit));
        jif3D::rvec InvCached1(CacheCalculator->LQDerivative(GravityTest, Misfit));
        jif3D::rvec InvCached2(CacheCalculator->LQDerivative(GravityTest, Misfit));
        BOOST_CHECK(InvDirect.size() == InvCached1.size());
        BOOST_CHECK(InvDirect.size() == InvCached2.size());
        for (size_t i = 0; i < InvDirect.size(); ++i)
          {
            BOOST_CHECK_CLOSE(InvDirect(i), InvCached1(i),
                std::numeric_limits<float>::epsilon());
            BOOST_CHECK_CLOSE(InvDirect(i), InvCached2(i),
                std::numeric_limits<float>::epsilon());
          }

      }

    BOOST_AUTO_TEST_CASE(diff_calc_test)
      {
        typedef typename jif3D::FullSensitivityGravMagCalculator<jif3D::ThreeDGravityModel> CalculatorType;
        boost::shared_ptr<CalculatorType> RawCalculator(
            jif3D::CreateGravityCalculator<CalculatorType>::MakeTensor());
        boost::shared_ptr<CalculatorType> InvCalculator(
            jif3D::CreateGravityCalculator<CalculatorType>::MakeTensor());
        boost::shared_ptr<jif3D::VectorTransform> Transform = boost::make_shared<
            jif3D::FTGInvariant>();
        InvCalculator->SetDataTransform(Transform);

        jif3D::ThreeDGravityModel GravityTest;
        const size_t nmeas = 10;
        const size_t ncells = 10;
        MakeRandomModel(GravityTest, ncells, nmeas, true);

        jif3D::rvec RawData(RawCalculator->Calculate(GravityTest));
        const jif3D::rmat RawSens(RawCalculator->GetSensitivities());

        jif3D::rvec InvData(InvCalculator->Calculate(GravityTest));

        const size_t nmod = RawSens.size2();
        jif3D::rmat InvarSens(nmeas, nmod);

        for (size_t i = 0; i < nmeas; ++i)
          {
            for (size_t j = 0; j < nmod; ++j)
              {
                InvarSens(i, j) = RawSens(i * 9, j) * RawData(i * 9 + 4)
                    + RawSens(i * 9 + 4, j) * RawData(i * 9)

                    + RawSens(i * 9 + 4, j) * RawData(i * 9 + 8)
                    + RawSens(i * 9 + 8, j) * RawData(i * 9 + 4)

                    + RawSens(i * 9, j) * RawData(i * 9 + 8)
                    + RawSens(i * 9 + 8, j) * RawData(i * 9)

                - 2.0 * RawSens(i * 9 + 3, j) * RawData(i * 9 + 3)
                    - 2.0 * RawSens(i * 9 + 7, j) * RawData(i * 9 + 7)
                    - 2.0 * RawSens(i * 9 + 2, j) * RawData(i * 9 + 2);
              }
          }

        jif3D::rvec Misfit(nmeas);
        std::fill(Misfit.begin(), Misfit.end(), 1.0);
        jif3D::rvec InvDeriv(InvCalculator->LQDerivative(GravityTest, Misfit));
        jif3D::rvec CompDeriv(2.0 * ublas::prod(trans(InvarSens), Misfit));
        for (size_t i = 0; i < InvDeriv.size(); ++i)
          {
            BOOST_CHECK_CLOSE(InvDeriv(i), CompDeriv(i),
                std::numeric_limits<float>::epsilon());
          }
      }

    BOOST_AUTO_TEST_CASE(finitediff_test)
      {
        jif3D::FTGInvariant InvarTrans;
        jif3D::rvec TestInput(9);
        std::generate(TestInput.begin(), TestInput.end(), rand);

        jif3D::rvec Invar(InvarTrans.Transform(TestInput));
        jif3D::rmat Deriv(InvarTrans.Derivative(TestInput));
        const double delta = 0.001;
        for (size_t i = 0; i < TestInput.size(); ++i)
          {
            jif3D::rvec DeltaInput(TestInput);
            DeltaInput(i) *= 1.0 + delta;
            double FD = (InvarTrans.Transform(DeltaInput)(0) - Invar(0))
                / (delta * TestInput(i));
            BOOST_CHECK_CLOSE(Deriv(0, i), FD, 0.01);
          }
      }
    BOOST_AUTO_TEST_SUITE_END()
