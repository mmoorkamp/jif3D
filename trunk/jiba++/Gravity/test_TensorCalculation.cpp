//============================================================================
// Name        : test_TensorCalculation.cpp
// Author      : Feb 11, 2009
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
#include "BasicGravElements.h"
#include "MinMemGravityCalculator.h"
#include "FullSensitivityGravityCalculator.h"
#include "ScalarOMPGravityImp.h"
#include "TensorOMPGravityImp.h"

BOOST_AUTO_TEST_SUITE( TensorGravity_Test_Suite )

//check some general properties of the FTG tensor that hold for any measurement above the surface
BOOST_AUTO_TEST_CASE  (random_tensor_test)
    {
      jiba::ThreeDGravityModel GravityTest;

      const size_t nmeas = 15;
      MakeRandomModel(GravityTest,nmeas);

      boost::shared_ptr<jiba::MinMemGravityCalculator> Calculator(jiba::CreateGravityCalculator<jiba::MinMemGravityCalculator>::MakeTensor());
      jiba::rvec tensormeas(Calculator->Calculate(GravityTest));
      for (size_t i = 0; i < nmeas; ++i)
        {
          // check that tensor is traceless
          BOOST_CHECK_CLOSE(tensormeas[i*9] + tensormeas[i*9+4],
              -tensormeas[i*9+8], std::numeric_limits<float>::epsilon());
          //check that tensor is symmetric
          BOOST_CHECK_CLOSE(tensormeas[i*9+1], tensormeas[i*9+3],
              std::numeric_limits<float>::epsilon());
          BOOST_CHECK_CLOSE(tensormeas[i*9+2], tensormeas[i*9+6],
              std::numeric_limits<float>::epsilon());
          BOOST_CHECK_CLOSE(tensormeas[i*9+5], tensormeas[i*9+7],
              std::numeric_limits<float>::epsilon());
        }
    }

  //compare the result from the tensor calculation with finite difference scalar calculations
  BOOST_AUTO_TEST_CASE(fd_tensor_test)
    {
      jiba::ThreeDGravityModel GravityTest;
      //we want to control the position of the measurements ourselves
      const size_t nmeas = 0;
      const double delta = 0.000001;
      //the value of delta in %
      const double precision = 0.05;
      MakeRandomModel(GravityTest, nmeas);
      //setup points for finite differencing in vertical and horizontal directions
      const double xpos = 50.0;
      const double ypos = 70.0;
      const double zpos = -5.1;
      GravityTest.AddMeasurementPoint(xpos, ypos, zpos);
      GravityTest.AddMeasurementPoint(xpos, ypos, zpos - delta);
      GravityTest.AddMeasurementPoint(xpos, ypos, zpos + delta);
      GravityTest.AddMeasurementPoint(xpos, ypos - delta, zpos);
      GravityTest.AddMeasurementPoint(xpos, ypos + delta, zpos);
      GravityTest.AddMeasurementPoint(xpos - delta, ypos, zpos);
      GravityTest.AddMeasurementPoint(xpos + delta, ypos, zpos);
      //perform the calculations
      boost::shared_ptr<jiba::MinMemGravityCalculator> ScalarCalculator(jiba::CreateGravityCalculator<jiba::MinMemGravityCalculator>::MakeScalar());
      boost::shared_ptr<jiba::MinMemGravityCalculator> TensorCalculator(jiba::CreateGravityCalculator<jiba::MinMemGravityCalculator>::MakeTensor());
      jiba::rvec tensormeas(
          TensorCalculator->Calculate(GravityTest));
      jiba::rvec scalarmeas(
          ScalarCalculator->Calculate(GravityTest));
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
      jiba::ThreeDGravityModel GravityTest;
      const size_t nmeas = 15;
      MakeRandomModel(GravityTest,nmeas);
      boost::shared_ptr<jiba::FullSensitivityGravityCalculator> TensorCalculator(jiba::CreateGravityCalculator<jiba::FullSensitivityGravityCalculator>::MakeTensor());

      //Calculate twice, once with normal calculation, once cached
      boost::posix_time::ptime startfirst =
      boost::posix_time::microsec_clock::local_time();
      jiba::rvec tensormeas1(
          TensorCalculator->Calculate(GravityTest));
      boost::posix_time::ptime endfirst =
      boost::posix_time::microsec_clock::local_time();
      jiba::rvec tensormeas2(
          TensorCalculator->Calculate(GravityTest));
      boost::posix_time::ptime endsecond =
      boost::posix_time::microsec_clock::local_time();
      //the second time it should be much faster
      BOOST_CHECK(endfirst - startfirst> (endsecond - endfirst));

      //check that invalidating the cache works
      TensorCalculator->SetSensitivities() *= 10;
      jiba::rvec meas3(TensorCalculator->Calculate(GravityTest));
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
      jiba::ThreeDGravityModel GravityTest;
      const size_t nhorcells = 10;
      const size_t nzcells = 10;
      GravityTest.SetXCellSizes().resize(boost::extents[nhorcells]);
      GravityTest.SetYCellSizes().resize(boost::extents[nhorcells]);
      GravityTest.SetZCellSizes().resize(boost::extents[nzcells]);
      for (size_t i = 0; i < nhorcells; ++i) // set the values of the inner cells

        {
          GravityTest.SetXCellSizes()[i] = 10;
          GravityTest.SetYCellSizes()[i] = 10;
        }
      std::fill_n(GravityTest.SetZCellSizes().origin(),nzcells,10.0);
      GravityTest.SetDensities().resize(
          boost::extents[nhorcells][nhorcells][nzcells]);
      const double density = 2.1;
      std::fill_n(GravityTest.SetDensities().origin(),GravityTest.SetDensities().num_elements(),density);

      GravityTest.AddMeasurementPoint(5,5, -1);
      GravityTest.AddMeasurementPoint(5,5,50);

      std::vector<double> bg_dens(1,density), bg_thick(1,100.0);
      GravityTest.SetBackgroundDensities(bg_dens);
      GravityTest.SetBackgroundThicknesses(bg_thick);
      boost::shared_ptr<jiba::MinMemGravityCalculator> TensorCalculator(jiba::CreateGravityCalculator<jiba::MinMemGravityCalculator>::MakeTensor());
      jiba::rvec TensorMeas(TensorCalculator->Calculate(GravityTest));
      const double PoissTerm = -4.0 * M_PI * jiba::Grav_const * density;

      const double Trace1 = TensorMeas(0) + TensorMeas(4) + TensorMeas(8);
      const double Trace2 = TensorMeas(9) + TensorMeas(13) + TensorMeas(17);

      BOOST_CHECK_CLOSE(PoissTerm, Trace1-Trace2,
          std::numeric_limits<float>::epsilon());
      BOOST_CHECK_CLOSE(density, Trace2/(4.0 * M_PI * jiba::Grav_const),
          std::numeric_limits<float>::epsilon());

    }

  BOOST_AUTO_TEST_CASE(tensor_cuda_test)
    {
      jiba::ThreeDGravityModel GravityTest;
      const size_t nmeas = 7;
      MakeRandomModel(GravityTest, nmeas);
      boost::shared_ptr<jiba::MinMemGravityCalculator> CPUCalculator(jiba::CreateGravityCalculator<jiba::MinMemGravityCalculator>::MakeTensor());
      boost::shared_ptr<jiba::MinMemGravityCalculator> CudaCalculator(jiba::CreateGravityCalculator<jiba::MinMemGravityCalculator>::MakeTensor(true));

      jiba::rvec cpumeas(
          CPUCalculator->Calculate(GravityTest));
      jiba::rvec cudameas(
          CudaCalculator->Calculate(GravityTest));
      BOOST_CHECK(cpumeas.size() == cudameas.size());
      for (size_t i = 0; i < cpumeas.size(); ++i)
        {
          BOOST_CHECK_CLOSE(cpumeas(i), cudameas(i),
              std::numeric_limits<float>::epsilon());

        }
    }

  BOOST_AUTO_TEST_CASE(tensor_cuda_caching_test)
    {
      jiba::ThreeDGravityModel GravityTest;
      const size_t nmeas = 20;
      MakeRandomModel(GravityTest, nmeas);
      boost::shared_ptr<jiba::FullSensitivityGravityCalculator> Calculator(jiba::CreateGravityCalculator<jiba::FullSensitivityGravityCalculator>::MakeTensor(true));

      jiba::rvec meas1(Calculator->Calculate(GravityTest));
      jiba::rvec meas2(Calculator->Calculate(GravityTest));

      for (size_t i = 0; i < meas1.size(); ++i)
        {
          BOOST_CHECK_CLOSE(meas1(i), meas2(i),
              std::numeric_limits<float>::epsilon());
        }
    }
  BOOST_AUTO_TEST_CASE(tensor_cuda_resize_test)
    {
      jiba::ThreeDGravityModel GravityTest;

      const size_t nmeas = 20;
      MakeRandomModel(GravityTest, nmeas);
      boost::shared_ptr<jiba::FullSensitivityGravityCalculator> Calculator(jiba::CreateGravityCalculator<jiba::FullSensitivityGravityCalculator>::MakeTensor(true));
      jiba::rvec meas1(Calculator->Calculate(GravityTest));
      GravityTest.ClearMeasurementPoints();
      MakeRandomModel(GravityTest,nmeas/2);
      BOOST_CHECK_NO_THROW( jiba::rvec meas2(Calculator->Calculate(GravityTest)));
    }

BOOST_AUTO_TEST_SUITE_END()
