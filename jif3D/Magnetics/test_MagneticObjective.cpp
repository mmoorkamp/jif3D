
//============================================================================
// Name        : test_GravityObjective.cpp
// Author      : Jun 6, 2009
// Version     :
// Copyright   : 2009, mmoorkamp
//============================================================================
#define BOOST_TEST_MODULE MagneticObjective test
#define BOOST_TEST_MAIN ...

#include "../Global/Jif3DTesting.h"
#include "../Gravity/test_common.h"
#include "../Inversion/ThreeDModelObjective.h"
#include "ThreeDSusceptibilityModel.h"
#include "ThreeDMagnetizationModel.h"
#include "TotalFieldMagneticData.h"
#include "ThreeComponentMagneticData.h"
#include "OMPMagneticSusceptibilityImp.h"
#include "OMPMagnetizationImp.h"
#include "MagneticTransforms.h"
#include "../GravMag/FullSensitivityGravMagCalculator.h"
#include "../GravMag/MinMemGravMagCalculator.h"
#include "../Global/Jif3DPlatformHelper.h"
#include <boost/test/included/unit_test.hpp>
#include <boost/test/tools/floating_point_comparison.hpp>
#include <boost/assign/std/vector.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/numeric/conversion/cast.hpp>
#include <boost/make_shared.hpp>



BOOST_AUTO_TEST_SUITE( MagneticObjective_Test_Suite )

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
            double FDGrad = (Objective.CalcMisfit(Forward)
                - Objective.CalcMisfit(Backward)) / (2 * delta);
            if (std::abs(FDGrad) > 1e-20)
              {
                BOOST_CHECK_CLOSE(FDGrad, Gradient(i), 0.01);
              }
            else
              {
                BOOST_CHECK(std::abs(Gradient(i)) < 1e-10);
                std::cout << FDGrad << " " << Gradient(i) << std::endl;
              }
          }
      }

    BOOST_AUTO_TEST_CASE (full_field_deriv_test)
      {
        jif3D::ThreeDSusceptibilityModel MagTest;
        jif3D::TotalFieldMagneticData Data;
        const size_t nmeas = 3;
        const size_t ncells = 5;
        MakeRandomModel(MagTest, Data, ncells, nmeas, false);
        double inclination = jif3D::platform::drand48();
        double declination = jif3D::platform::drand48();
        double fieldstrength = 1.0 + jif3D::platform::drand48();

        typedef typename jif3D::MinMemGravMagCalculator<jif3D::TotalFieldMagneticData> CalculatorType;
        boost::shared_ptr<jif3D::ThreeDGravMagImplementation<jif3D::TotalFieldMagneticData> > Implementation(
            new jif3D::OMPMagneticSusceptibilityImp(inclination, declination, fieldstrength));

        boost::shared_ptr<CalculatorType> Calculator(new CalculatorType(Implementation));
        jif3D::rvec Observed(Calculator->Calculate(MagTest, Data));

        Observed *= 1.1;
        std::vector<double> d(Observed.size()), e(Observed.size(), 1.0);
        std::copy(Observed.begin(), Observed.end(), d.begin());
        Data.SetDataAndErrors(d, e);

        jif3D::ThreeDModelObjective<CalculatorType> Objective(*Calculator.get());
        Objective.SetObservedData(Data);
        Objective.SetCoarseModelGeometry(MagTest);
        jif3D::rvec InvModel(MagTest.GetSusceptibilities().num_elements());
        std::copy(MagTest.GetSusceptibilities().origin(),
            MagTest.GetSusceptibilities().origin()
                + MagTest.GetSusceptibilities().num_elements(), InvModel.begin());

        CheckGradient(Objective, InvModel);
      }

    BOOST_AUTO_TEST_CASE (T_deriv_test)
      {
        jif3D::ThreeDSusceptibilityModel MagTest;
        jif3D::TotalFieldMagneticData Data;

        const size_t nmeas = 3;
        const size_t ncells = 5;
        MakeRandomModel(MagTest, Data, ncells, nmeas, false);
        double inclination = jif3D::platform::drand48();
        double declination = jif3D::platform::drand48();
        double fieldstrength = 1.0 + jif3D::platform::drand48();

        typedef typename jif3D::MinMemGravMagCalculator<jif3D::TotalFieldMagneticData> CalculatorType;
        boost::shared_ptr<jif3D::ThreeDGravMagImplementation<jif3D::TotalFieldMagneticData> > Implementation(
            new jif3D::OMPMagneticSusceptibilityImp(inclination, declination, fieldstrength));

        boost::shared_ptr<CalculatorType> Calculator(new CalculatorType(Implementation));
        Calculator->SetDataTransform(boost::make_shared<jif3D::TotalField>());
        jif3D::rvec Observed(Calculator->Calculate(MagTest, Data));

        Observed *= 1.1;
        std::vector<double> d(Observed.size()), e(Observed.size(), 1.0);
        std::copy(Observed.begin(), Observed.end(), d.begin());
        Data.SetDataAndErrors(d, e);

        jif3D::ThreeDModelObjective<CalculatorType> Objective(*Calculator.get());
        Objective.SetDataTransform(boost::make_shared<jif3D::TotalField>());
        Objective.SetObservedData(Data);
        Objective.SetCoarseModelGeometry(MagTest);
        jif3D::rvec InvModel(MagTest.GetSusceptibilities().num_elements());
        std::copy(MagTest.GetSusceptibilities().origin(),
            MagTest.GetSusceptibilities().origin()
                + MagTest.GetSusceptibilities().num_elements(), InvModel.begin());

        CheckGradient(Objective, InvModel);
      }

    BOOST_AUTO_TEST_CASE (deltaT_deriv_test)
      {
        jif3D::ThreeDSusceptibilityModel MagTest;
        jif3D::TotalFieldMagneticData Data;

        const size_t nmeas = 3;
        const size_t ncells = 5;
        MakeRandomModel(MagTest, Data, ncells, nmeas, false);
        double inclination = jif3D::platform::drand48();
        double declination = jif3D::platform::drand48();
        double fieldstrength = 1.0 + jif3D::platform::drand48();

        typedef typename jif3D::MinMemGravMagCalculator<jif3D::TotalFieldMagneticData> CalculatorType;
        boost::shared_ptr<jif3D::ThreeDGravMagImplementation<jif3D::TotalFieldMagneticData> > Implementation(
            new jif3D::OMPMagneticSusceptibilityImp(inclination, declination, fieldstrength));

        boost::shared_ptr<CalculatorType> Calculator(new CalculatorType(Implementation));
        Calculator->SetDataTransform(
            boost::make_shared<jif3D::TotalFieldAnomaly>(inclination, declination));

        jif3D::rvec Observed(Calculator->Calculate(MagTest, Data));

        Observed *= 1.1;
        std::vector<double> d(Observed.size()), e(Observed.size(), 1.0);
        std::copy(Observed.begin(), Observed.end(), d.begin());
        Data.SetDataAndErrors(d, e);

        jif3D::ThreeDModelObjective<CalculatorType> Objective(*Calculator.get());
        Objective.SetDataTransform(
            boost::make_shared<jif3D::TotalFieldAnomaly>(inclination, declination));
        Objective.SetObservedData(Data);
        Objective.SetCoarseModelGeometry(MagTest);
        jif3D::rvec InvModel(MagTest.GetSusceptibilities().num_elements());
        std::copy(MagTest.GetSusceptibilities().origin(),
            MagTest.GetSusceptibilities().origin()
                + MagTest.GetSusceptibilities().num_elements(), InvModel.begin());

        CheckGradient(Objective, InvModel);
      }


    BOOST_AUTO_TEST_CASE (vector_deriv_test)
      {
        jif3D::ThreeDMagnetizationModel MagTest;
        jif3D::ThreeComponentMagneticData Data;
        const size_t nmeas = 3;
        const size_t ncells = 5;
        MakeRandomModel(MagTest, Data, ncells, nmeas, false);
        MagTest.SetMagnetization_Y() = MagTest.GetMagnetization_X();
        MagTest.SetMagnetization_Z() = MagTest.GetMagnetization_X();


        typedef typename jif3D::MinMemGravMagCalculator<jif3D::ThreeComponentMagneticData> CalculatorType;
        boost::shared_ptr<jif3D::ThreeDGravMagImplementation<jif3D::ThreeComponentMagneticData> > Implementation(
            new jif3D::OMPMagnetizationImp());

        boost::shared_ptr<CalculatorType> Calculator(new CalculatorType(Implementation));
        jif3D::rvec Observed(Calculator->Calculate(MagTest, Data));
        std::cout << Observed << std::endl;
        Observed *= 1.1;
        std::vector<double> d(Observed.size()), e(Observed.size(), 1.0);
        std::copy(Observed.begin(), Observed.end(), d.begin());
        Data.SetDataAndErrors(d, e);

        jif3D::ThreeDModelObjective<CalculatorType> Objective(*Calculator.get());
        Objective.SetObservedData(Data);
        Objective.SetCoarseModelGeometry(MagTest);
        jif3D::rvec InvModel(MagTest.GetModelParameters());
        std::cout << MagTest.GetModelShape() << std::endl;

        std::cout << InvModel << std::endl;

        CheckGradient(Objective, InvModel);
      }

    BOOST_AUTO_TEST_SUITE_END()
