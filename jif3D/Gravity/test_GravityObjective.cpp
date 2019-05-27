//============================================================================
// Name        : test_GravityObjective.cpp
// Author      : Jun 6, 2009
// Version     :
// Copyright   : 2009, mmoorkamp
//============================================================================

#define BOOST_TEST_MODULE ThreeDGravityModel test
#define BOOST_TEST_MAIN ...
#include "../Global/Jif3DTesting.h"
#include <boost/test/included/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/assign/std/vector.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/numeric/conversion/cast.hpp>

#include "test_common.h"
#include "../Inversion/ThreeDModelObjective.h"
#include "ThreeDGravityFactory.h"
#include "ThreeDGravityModel.h"
#include "ScalarGravityData.h"
#include "TensorGravityData.h"
#include "../GravMag/FullSensitivityGravMagCalculator.h"
#include "../GravMag/MinMemGravMagCalculator.h"

BOOST_AUTO_TEST_SUITE (GravityObjective_Test_Suite)

//we check the gradient from the forward modeling by comparing it to the result
//of a finite difference computation, this is the same for FTG and scalar data
//and therefore in a separate function
void CheckGradient(jif3D::ObjectiveFunction &Objective, const jif3D::rvec &Model)
  {
    Objective.CalcMisfit(Model);
    jif3D::rvec Gradient = Objective.CalcGradient(Model);

    BOOST_CHECK(Model.size() == Gradient.size());
    for (size_t i = 0; i < Gradient.size(); ++i)
      {
        double delta = Model(i) * 0.01;
        jif3D::rvec Forward(Model);
        jif3D::rvec Backward(Model);
        Forward(i) += delta;
        Backward(i) -= delta;
        if (delta > 0.0)
          {
            double FDGrad = (Objective.CalcMisfit(Forward)
                - Objective.CalcMisfit(Backward)) / (2 * delta);
            BOOST_CHECK_CLOSE(FDGrad, Gradient(i), 0.001);
          }
      }
  }

//check the derivative for the FTG tensor
BOOST_AUTO_TEST_CASE (tensor_fdderiv_test)
  {
    jif3D::ThreeDGravityModel GravityTest;
    jif3D::TensorGravityData Data;
    const size_t nmeas = 3;
    const size_t ncells = 5;

    MakeRandomModel(GravityTest,Data, ncells,nmeas,false);
    typedef typename jif3D::FullSensitivityGravMagCalculator<jif3D::TensorGravityData> CalculatorType;
    boost::shared_ptr<CalculatorType> TensorCalculator(
        jif3D::CreateGravityCalculator<CalculatorType>::MakeTensor());
    jif3D::rvec Observed(TensorCalculator->Calculate(GravityTest, Data));
    //we have to have different data from our observed data
    //otherwise our gradient will be zero
    Observed *= 1.1;
    std::vector<double> d(Observed.size()), e(Observed.size(), 1.0);
    std::copy(Observed.begin(),Observed.end(),d.begin());
    Data.SetDataAndErrors(d,e);
    jif3D::ThreeDModelObjective<CalculatorType> FTGObjective(*TensorCalculator.get());

    FTGObjective.SetObservedData(Data);
    FTGObjective.SetCoarseModelGeometry(GravityTest);
    jif3D::rvec InvModel(GravityTest.GetDensities().num_elements());
    std::copy(GravityTest.GetDensities().origin(),
        GravityTest.GetDensities().origin()
        + GravityTest.GetDensities().num_elements(), InvModel.begin());

    CheckGradient(FTGObjective,InvModel);
  }

BOOST_AUTO_TEST_CASE (scalar_fdderiv_test)
  {
    jif3D::ThreeDGravityModel GravityTest;
    jif3D::ScalarGravityData Data;
    const size_t nmeas = 3;
    const size_t ncells = 5;
    MakeRandomModel(GravityTest,Data, ncells,nmeas,false);
    typedef typename jif3D::FullSensitivityGravMagCalculator<jif3D::ScalarGravityData> CalculatorType;
    boost::shared_ptr<CalculatorType> ScalarCalculator(
        jif3D::CreateGravityCalculator<CalculatorType>::MakeScalar());
    jif3D::rvec Observed(ScalarCalculator->Calculate(GravityTest, Data));
    Observed *= 1.1;
    std::vector<double> d(Observed.size()), e(Observed.size(), 1.0);
    std::copy(Observed.begin(),Observed.end(),d.begin());
    Data.SetDataAndErrors(d,e);


    jif3D::ThreeDModelObjective<CalculatorType> ScalarObjective(*ScalarCalculator.get());
    ScalarObjective.SetObservedData(Data);
    ScalarObjective.SetCoarseModelGeometry(GravityTest);
    jif3D::rvec InvModel(GravityTest.GetDensities().num_elements());
    std::copy(GravityTest.GetDensities().origin(),
        GravityTest.GetDensities().origin()
        + GravityTest.GetDensities().num_elements(), InvModel.begin());

    CheckGradient(ScalarObjective,InvModel);
  }

BOOST_AUTO_TEST_SUITE_END()
