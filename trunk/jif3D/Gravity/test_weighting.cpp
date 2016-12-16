//============================================================================
// Name        : test_weighting.cpp
// Author      : Sep 18, 2008
// Version     :
// Copyright   : 2008, mmoorkamp
//============================================================================


#define BOOST_TEST_MODULE GravityDepthWeighting test
#define BOOST_TEST_MAIN ...
#include "../Global/Jif3DTesting.h"
#include <boost/test/included/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <stdlib.h>
#include <time.h>
#include "DepthWeighting.h"
#include "ThreeDGravityFactory.h"
#include "../GravMag/ThreeDGravMagCalculator.h"
#include "../GravMag/FullSensitivityGravMagCalculator.h"

BOOST_AUTO_TEST_SUITE( Gravity_DepthWeighting_Suite )

BOOST_AUTO_TEST_CASE(weightingderiv_test)
    {
      //this value is arbitrary
      const double z0 = acos(-1.0);
      const double exponent = -2.0;
      const double z = 0.4;
      //we compare the analytical derivative to a numerical derivative
      double value1 = jif3D::WeightingTerm(exponent)(z, z0);
      const double delta = std::numeric_limits<float>::epsilon();
      double value2 = jif3D::WeightingTerm(exponent)(z + 2.0 * delta, z0);
      double deriv = jif3D::WeightingTerm(exponent).deriv(z + delta, z0);
      BOOST_CHECK_CLOSE(deriv, (value2-value1)/(2.0*delta), 0.1);
    }
  BOOST_AUTO_TEST_CASE(weightingaverage_test)
    {
      //we want to check the correctness of the average
      //we choose z0, and the two boundaries for the average
      //so that the function is relatively constant over the intervall
      //that makes it easy to compare to a numerical average
      const double z0 = 1.0;
      const double exponent = -3.0;
      const double z1 = 200.0;
      const double z2 = 201.0;
      //out analytical average
      const double WeightingAverage = jif3D::WeightingTerm(exponent).average(z1,
          z2, z0);
      const double eps = 0.001;
      const int steps = (z2 - z1) / eps;
      double myavg = 0.0;
      //sum up values in the interval and divide by the number of interval
      for (int i = 0; i < steps; ++i)
        myavg += jif3D::WeightingTerm(exponent)(z1 + i * eps, z0);
      myavg /= steps;
      BOOST_CHECK_CLOSE(WeightingAverage, myavg, 0.001);
    }
  BOOST_AUTO_TEST_CASE(fitz0_test)
    {
      jif3D::ThreeDGravityModel Model;
      const size_t nz = 20;
      //we only need the model to make FitZ0 happy
      //all it does is determine the starting index
      //we are only interested in z dependence
      //so we have only 1 cell in x and y
      Model.SetXCellSizes().resize(boost::extents[1]);
      Model.SetYCellSizes().resize(boost::extents[1]);
      Model.SetZCellSizes().resize(boost::extents[nz]);
      Model.SetXCellSizes()[0] = 10.0;
      Model.SetYCellSizes()[0] = 10.0;
      Model.AddMeasurementPoint(5, 5, -1);
      //we make up some "sensitivities"
      jif3D::rvec PseudoSens(nz);
      //we want to get back this number
      const double z0 = 237.0;
      //this is the exponent we use to fit the gravity data
      const double exponent = -3.0;
      double currdepth = 0.0;
      //fill cell sizes and "sensitivities" with well defined values
      for (size_t i = 0; i < nz; ++i)
        {
          Model.SetZCellSizes()[i] = nz * i;
          currdepth += Model.GetZCellSizes()[i];
          PseudoSens[i] = jif3D::WeightingTerm(exponent)(currdepth, z0);
        }
      //FitZ0 should find z0 specified above
      double foundz0 = jif3D::FitZ0(PseudoSens, Model.GetZCellSizes(),
          jif3D::WeightingTerm(exponent));
      //in reality the gradient is small at the end, so we can have 1% difference
      BOOST_CHECK_CLOSE(z0, foundz0, 1);
    }

  BOOST_AUTO_TEST_CASE(extract_test)
    {
      jif3D::ThreeDGravityModel Model;

      const size_t ncells = 3;
      Model.SetXCellSizes().resize(boost::extents[ncells]);
      Model.SetYCellSizes().resize(boost::extents[ncells]);
      Model.SetZCellSizes().resize(boost::extents[ncells]);
      for (size_t i = 0; i < ncells; ++i)
        {
          Model.SetXCellSizes()[i] = 10.0;
          Model.SetYCellSizes()[i] = 10.0;
          Model.SetZCellSizes()[i] = 10.0;
        }
      Model.AddMeasurementPoint(2, 2, -1);
      Model.AddMeasurementPoint(15, 15, -1);
      Model.AddMeasurementPoint(28, 28, -1);
      Model.SetDensities().resize(boost::extents[ncells][ncells][ncells]);
      const double density = 2.1;
      std::fill_n(Model.SetDensities().origin(),
          Model.SetDensities().num_elements(), density);
      boost::shared_ptr<jif3D::FullSensitivityGravMagCalculator>
          Calculator(jif3D::CreateGravityCalculator<
              jif3D::FullSensitivityGravMagCalculator>::MakeTensor());
      Calculator->Calculate(Model);
      jif3D::rvec SensProfile;
      jif3D::ExtractMiddleSens(Model, Calculator->GetSensitivities(),
          Calculator->GetDataPerMeasurement(), SensProfile);
      size_t offset = Model.IndexToOffset(1, 1, 0);
      jif3D::rvec Compare(ncells);
      size_t startindex = offset;
      boost::numeric::ublas::matrix_row<const jif3D::rmat> Row(Calculator->GetSensitivities(),Calculator->GetDataPerMeasurement());
      std::copy(Row.begin() + startindex,
          Row.begin() + startindex + ncells,
          Compare.begin());
      for (size_t i = 0; i < ncells; ++i)
        {
          BOOST_CHECK_CLOSE(Compare(i),SensProfile(i),std::numeric_limits<float>::epsilon());
        }
      BOOST_CHECK(std::equal(Compare.begin(),Compare.end(),SensProfile.begin()));
    }
BOOST_AUTO_TEST_SUITE_END()
