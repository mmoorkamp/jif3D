//============================================================================
// Name        : test_TensorCalculation.cpp
// Author      : Feb 11, 2009
// Version     :
// Copyright   : 2009, mmoorkamp
//============================================================================

#define BOOST_TEST_MODULE ThreeDGravityModel test
#define BOOST_TEST_MAIN ...
#include "../Global/Jif3DTesting.h"
#include <boost/test/included/unit_test.hpp>
#include <boost/test/tools/floating_point_comparison.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/numeric/conversion/cast.hpp>
#include <boost/math/constants/constants.hpp>

#include "test_common.h"
#include "../Global/Jif3DPlatformHelper.h"
#include "BasicGravElements.h"
#include "../GravMag/MinMemGravMagCalculator.h"
#include "../GravMag/FullSensitivityGravMagCalculator.h"
#include "../GravMag/DiskGravMagCalculator.h"
#include "ScalarOMPGravityImp.h"
#include "TensorOMPGravityImp.h"
#include "ThreeDGravityFactory.h"
#include "TensorGravityData.h"
#include "ScalarGravityData.h"

BOOST_AUTO_TEST_SUITE (TensorGravity_Test_Suite)

//check some general properties of the FTG tensor that hold for any measurement above the surface
    BOOST_AUTO_TEST_CASE (random_tensor_test)
      {
        jif3D::ThreeDGravityModel GravityTest;
        jif3D::TensorGravityData Data;
        const size_t nmeas = 15;
        const size_t ncells = 15;
        MakeRandomModel(GravityTest, Data, ncells, nmeas);

        typedef typename jif3D::FullSensitivityGravMagCalculator<jif3D::TensorGravityData> CalculatorType;
        boost::shared_ptr<CalculatorType> Calculator(
            jif3D::CreateGravityCalculator<CalculatorType>::MakeTensor());
        jif3D::rvec tensormeas(Calculator->Calculate(GravityTest, Data));
        for (size_t i = 0; i < nmeas; ++i)
          {
            // check that tensor is traceless
            BOOST_CHECK_CLOSE(tensormeas[i * 9] + tensormeas[i * 9 + 4],
                -tensormeas[i * 9 + 8], std::numeric_limits<float>::epsilon());
            //check that tensor is symmetric
            BOOST_CHECK_CLOSE(tensormeas[i * 9 + 1], tensormeas[i * 9 + 3],
                std::numeric_limits<float>::epsilon());
            BOOST_CHECK_CLOSE(tensormeas[i * 9 + 2], tensormeas[i * 9 + 6],
                std::numeric_limits<float>::epsilon());
            BOOST_CHECK_CLOSE(tensormeas[i * 9 + 5], tensormeas[i * 9 + 7],
                std::numeric_limits<float>::epsilon());
          }
      }

//compare the result from the tensor calculation with finite difference scalar calculations
    BOOST_AUTO_TEST_CASE(fd_tensor_test)
      {
        jif3D::ThreeDGravityModel GravityTest;
        jif3D::TensorGravityData TensorData;
        jif3D::ScalarGravityData ScalarData;

        //we want to control the position of the measurements ourselves
        const size_t nmeas = 0;
        const size_t ncells = 10;
        const double delta = 0.000001;
        //the value of delta in %
        const double precision = 0.05;
        MakeRandomModel(GravityTest, TensorData, ncells, nmeas);
        //setup points for finite differencing in vertical and horizontal directions
        const double xpos = 50.0;
        const double ypos = 70.0;
        const double zpos = -5.1;
        TensorData.AddMeasurementPoint(xpos, ypos, zpos);
        TensorData.AddMeasurementPoint(xpos, ypos, zpos - delta);
        TensorData.AddMeasurementPoint(xpos, ypos, zpos + delta);
        TensorData.AddMeasurementPoint(xpos, ypos - delta, zpos);
        TensorData.AddMeasurementPoint(xpos, ypos + delta, zpos);
        TensorData.AddMeasurementPoint(xpos - delta, ypos, zpos);
        TensorData.AddMeasurementPoint(xpos + delta, ypos, zpos);
        ScalarData.CopyMeasurementConfiguration(TensorData);
        //perform the calculations
        typedef typename jif3D::MinMemGravMagCalculator<jif3D::ScalarGravityData> ScalarCalculatorType;
        typedef typename jif3D::MinMemGravMagCalculator<jif3D::TensorGravityData> TensorCalculatorType;
        boost::shared_ptr<ScalarCalculatorType> ScalarCalculator(
            jif3D::CreateGravityCalculator<ScalarCalculatorType>::MakeScalar());
        boost::shared_ptr<TensorCalculatorType> TensorCalculator(
            jif3D::CreateGravityCalculator<TensorCalculatorType>::MakeTensor());
        jif3D::rvec tensormeas(TensorCalculator->Calculate(GravityTest, TensorData));
        jif3D::rvec scalarmeas(ScalarCalculator->Calculate(GravityTest, ScalarData));
        //extract the tensor components
        double uzz = tensormeas(8);
        double uzy = tensormeas(7);
        double uzx = tensormeas(6);
        //calculate the same components through finite differencing
        double fduzz = (scalarmeas(2) - scalarmeas(1)) / (2 * delta);
        double fduzy = (scalarmeas(4) - scalarmeas(3)) / (2 * delta);
        double fduzx = (scalarmeas(6) - scalarmeas(5)) / (2 * delta);
        //check that they match within the precision of delta
        BOOST_CHECK_CLOSE(uzz, fduzz, precision);
        BOOST_CHECK_CLOSE(uzy, fduzy, precision);
        BOOST_CHECK_CLOSE(uzx, fduzx, precision);
      }

//check whether calculation by the sensitivity matrix yields the same results
//for tensor calculations
    BOOST_AUTO_TEST_CASE(tensor_caching_test)
      {
        jif3D::ThreeDGravityModel GravityTest;
        jif3D::TensorGravityData TensorData;

        const size_t nmeas = 15;
        MakeRandomModel(GravityTest, TensorData, nmeas);
        typedef typename jif3D::FullSensitivityGravMagCalculator<jif3D::TensorGravityData> CalculatorType;
        boost::shared_ptr<CalculatorType> TensorCalculator(
            jif3D::CreateGravityCalculator<CalculatorType>::MakeTensor());

        //Calculate twice, once with normal calculation, once cached
        boost::posix_time::ptime startfirst =
            boost::posix_time::microsec_clock::local_time();
        jif3D::rvec tensormeas1(TensorCalculator->Calculate(GravityTest, TensorData));
        boost::posix_time::ptime endfirst =
            boost::posix_time::microsec_clock::local_time();
        jif3D::rvec tensormeas2(TensorCalculator->Calculate(GravityTest, TensorData));
        boost::posix_time::ptime endsecond =
            boost::posix_time::microsec_clock::local_time();
        //the second time it should be much faster
        BOOST_CHECK(endfirst - startfirst > (endsecond - endfirst));

        //check that invalidating the cache works
        TensorCalculator->SetSensitivities() *= 10;
        jif3D::rvec meas3(TensorCalculator->Calculate(GravityTest, TensorData));
        //compare all tensor elements
        for (size_t i = 0; i < tensormeas1.size(); ++i)
          {
            BOOST_CHECK_CLOSE(tensormeas1(i), tensormeas2(i),
                std::numeric_limits<float>::epsilon());
            BOOST_CHECK_CLOSE(tensormeas1(i), meas3(i),
                std::numeric_limits<float>::epsilon());

          }
      }

//check whether calculation by the sensitivity matrix yields the same results
//for tensor calculations
    BOOST_AUTO_TEST_CASE(tensor_disk_test)
      {
        jif3D::ThreeDGravityModel GravityTest;
        jif3D::TensorGravityData TensorData;

        const size_t nmeas = 5;
        MakeRandomModel(GravityTest, TensorData, nmeas);
        typedef typename jif3D::DiskGravMagCalculator<jif3D::TensorGravityData> CalculatorType;
        boost::shared_ptr<CalculatorType> TensorCalculator(
            jif3D::CreateGravityCalculator<CalculatorType>::MakeTensor());

        //Calculate twice, once with normal calculation, once cached
        boost::posix_time::ptime startfirst =
            boost::posix_time::microsec_clock::local_time();
        jif3D::rvec tensormeas1(TensorCalculator->Calculate(GravityTest, TensorData));
        boost::posix_time::ptime endfirst =
            boost::posix_time::microsec_clock::local_time();
        jif3D::rvec tensormeas2(TensorCalculator->Calculate(GravityTest, TensorData));
        boost::posix_time::ptime endsecond =
            boost::posix_time::microsec_clock::local_time();
        //the second time it should be much faster
        BOOST_CHECK(endfirst - startfirst > (endsecond - endfirst));

        //check that invalidating the cache works
        TensorCalculator->SetSensitivities() *= 10;
        jif3D::rvec meas3(TensorCalculator->Calculate(GravityTest, TensorData));
        //compare all tensor elements
        for (size_t i = 0; i < tensormeas1.size(); ++i)
          {
            BOOST_CHECK_CLOSE(tensormeas1(i), tensormeas2(i),
                std::numeric_limits<float>::epsilon());
            BOOST_CHECK_CLOSE(tensormeas1(i), meas3(i),
                std::numeric_limits<float>::epsilon());

          }
      }

//check whether the diagonal elements of the tensor obey the poisson equation
    BOOST_AUTO_TEST_CASE(tensor_poisson_test)
      {
        jif3D::ThreeDGravityModel GravityTest;
        jif3D::TensorGravityData TensorData;

        const size_t nhorcells = 10;
        const size_t nzcells = 10;
        GravityTest.SetMeshSize(nhorcells, nhorcells, nzcells);
        jif3D::ThreeDModelBase::t3DModelDim XCD(nhorcells, 10.0), YCD(nhorcells, 10.0),
            ZCD(nzcells, 10.0);
        GravityTest.SetXCellSizes(XCD);
        GravityTest.SetYCellSizes(YCD);
        GravityTest.SetZCellSizes(ZCD);

        const double density = 2.1;
        std::fill_n(GravityTest.SetDensities().origin(),
            GravityTest.SetDensities().num_elements(), density);

        TensorData.AddMeasurementPoint(5, 5, -1);
        TensorData.AddMeasurementPoint(5, 5, 50);

        std::vector<double> bg_dens(1, density), bg_thick(1, 100.0);
        GravityTest.SetBackgroundDensities(bg_dens);
        GravityTest.SetBackgroundThicknesses(bg_thick);
        typedef typename jif3D::MinMemGravMagCalculator<jif3D::TensorGravityData> CalculatorType;
        boost::shared_ptr<CalculatorType> TensorCalculator(
            jif3D::CreateGravityCalculator<CalculatorType>::MakeTensor());
        jif3D::rvec TensorMeas(TensorCalculator->Calculate(GravityTest, TensorData));
        const double PoissTerm = -4.0 * boost::math::constants::pi<double>()
            * jif3D::Grav_const * density;

        const double Trace1 = TensorMeas(0) + TensorMeas(4) + TensorMeas(8);
        const double Trace2 = TensorMeas(9) + TensorMeas(13) + TensorMeas(17);

        BOOST_CHECK_CLOSE(PoissTerm, Trace1 - Trace2,
            std::numeric_limits<float>::epsilon());
        BOOST_CHECK_CLOSE(density,
            Trace2 / (4.0 * boost::math::constants::pi<double>() * jif3D::Grav_const),
            std::numeric_limits<float>::epsilon());

      }
#ifdef HAVEGPU
    BOOST_AUTO_TEST_CASE(tensor_cuda_test)
      {
        std::cout << " Running CUDA test " << std::endl;
        jif3D::ThreeDGravityModel GravityTest;
        jif3D::TensorGravityData TensorData;

        const size_t nmeas = 7;
        const size_t ncells = 7;
        MakeRandomModel(GravityTest,TensorData, ncells, nmeas);
        typedef typename jif3D::MinMemGravMagCalculator<jif3D::TensorGravityData> CalculatorType;
        boost::shared_ptr<CalculatorType> CPUCalculator(jif3D::CreateGravityCalculator<CalculatorType>::MakeTensor());
        boost::shared_ptr<CalculatorType> CudaCalculator(jif3D::CreateGravityCalculator<CalculatorType>::MakeTensor(true));

        jif3D::rvec cpumeas(
            CPUCalculator->Calculate(GravityTest,TensorData));
        jif3D::rvec cudameas(
            CudaCalculator->Calculate(GravityTest,TensorData));
        BOOST_CHECK(cpumeas.size() == cudameas.size());
        for (size_t i = 0; i < cpumeas.size(); ++i)
          {
            BOOST_CHECK_CLOSE(cpumeas(i), cudameas(i),
                std::numeric_limits<float>::epsilon());

          }
        const size_t nruns = 10;
        for (size_t i = 0; i < nruns; ++i)
          {
            jif3D::rvec newcuda(CudaCalculator->Calculate(GravityTest,TensorData));
            BOOST_CHECK(std::equal(newcuda.begin(),newcuda.end(),cudameas.begin()));
          }
      }

    BOOST_AUTO_TEST_CASE(tensor_cuda_caching_test)
      {
        jif3D::ThreeDGravityModel GravityTest;
        jif3D::TensorGravityData TensorData;

        const size_t nmeas = 20;
        const size_t ncells = 7;

        MakeRandomModel(GravityTest, TensorData, ncells, nmeas);
        typedef typename jif3D::FullSensitivityGravMagCalculator<jif3D::TensorGravityData> CalculatorType;
        boost::shared_ptr<CalculatorType> Calculator(jif3D::CreateGravityCalculator<CalculatorType>::MakeTensor(true));

        jif3D::rvec meas1(Calculator->Calculate(GravityTest,TensorData));
        jif3D::rvec meas2(Calculator->Calculate(GravityTest,TensorData));

        for (size_t i = 0; i < meas1.size(); ++i)
          {
            BOOST_CHECK_CLOSE(meas1(i), meas2(i),
                std::numeric_limits<float>::epsilon());
          }
      }

#endif

    BOOST_AUTO_TEST_CASE (lqderivative_test)
      {
        jif3D::ThreeDGravityModel GravityTest;
        jif3D::TensorGravityData TensorData;
        const size_t ncells = 10;
        const size_t nmeas = 10;
        MakeRandomModel(GravityTest, TensorData, ncells, nmeas, false);
        typedef typename jif3D::FullSensitivityGravMagCalculator<jif3D::TensorGravityData> CalculatorType;
        boost::shared_ptr<CalculatorType> TensorCalculator(
            jif3D::CreateGravityCalculator<CalculatorType>::MakeTensor());
        jif3D::rvec Misfit(nmeas * TensorCalculator->GetDataPerMeasurement());
        std::generate(Misfit.begin(), Misfit.end(), jif3D::platform::drand48);
        jif3D::rvec Deriv(
            TensorCalculator->LQDerivative(GravityTest, TensorData, Misfit));
        TensorCalculator->Calculate(GravityTest, TensorData);
        jif3D::rvec Compare(
            2.0
                * boost::numeric::ublas::prec_prod(
                    ublas::trans(TensorCalculator->GetSensitivities()), Misfit));
        //and test the caching, too
        jif3D::rvec Deriv2(
            TensorCalculator->LQDerivative(GravityTest, TensorData, Misfit));
        const size_t ngrid = GravityTest.GetDensities().num_elements();
        BOOST_CHECK(Deriv.size() == ngrid);
        BOOST_CHECK(Compare.size() == ngrid);
        BOOST_CHECK(Deriv2.size() == ngrid);
        for (size_t i = 0; i < ngrid; ++i)
          {
            BOOST_CHECK_CLOSE(Deriv(i), Compare(i),
                std::numeric_limits<float>::epsilon());
            BOOST_CHECK_CLOSE(Deriv2(i), Compare(i),
                std::numeric_limits<float>::epsilon());
          }
      }

    BOOST_AUTO_TEST_SUITE_END()
