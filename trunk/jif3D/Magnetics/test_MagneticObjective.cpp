//============================================================================
// Name        : test_GravityObjective.cpp
// Author      : Jun 6, 2009
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

#include "../Gravity/test_common.h"
#include "../Inversion/ThreeDModelObjective.h"
#include "ThreeDMagneticModel.h"
#include "OMPMagneticImp.h"
#include "MagneticTransforms.h"
#include "../GravMag/FullSensitivityGravMagCalculator.h"
#include "../GravMag/MinMemGravMagCalculator.h"

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
            if (delta > 0.0)
              {
                double FDGrad = (Objective.CalcMisfit(Forward)
                    - Objective.CalcMisfit(Backward)) / (2 * delta);
                BOOST_CHECK_CLOSE(FDGrad, Gradient(i), 0.001);
              }
          }
      }

    BOOST_AUTO_TEST_CASE (full_field_deriv_test)
      {
        jif3D::ThreeDMagneticModel MagTest;
        const size_t nmeas = 3;
        const size_t ncells = 5;
        MakeRandomModel(MagTest, ncells, nmeas, false);
        double inclination = drand48();
        double declination = drand48();
        double fieldstrength = 1.0 + drand48();

        typedef typename jif3D::MinMemGravMagCalculator<jif3D::ThreeDMagneticModel> CalculatorType;
        boost::shared_ptr<jif3D::ThreeDGravMagImplementation<jif3D::ThreeDMagneticModel> > Implementation(
            new jif3D::OMPMagneticImp(inclination, declination, fieldstrength));

        boost::shared_ptr<CalculatorType> Calculator(new CalculatorType(Implementation));
        jif3D::rvec Observed(Calculator->Calculate(MagTest));

        Observed *= 1.1;

        jif3D::ThreeDModelObjective<CalculatorType> Objective(*Calculator.get());
        Objective.SetObservedData(Observed);
        Objective.SetCoarseModelGeometry(MagTest);
        jif3D::rvec InvModel(MagTest.GetSusceptibilities().num_elements());
        std::copy(MagTest.GetSusceptibilities().origin(),
            MagTest.GetSusceptibilities().origin()
                + MagTest.GetSusceptibilities().num_elements(), InvModel.begin());

        CheckGradient(Objective, InvModel);
      }

    BOOST_AUTO_TEST_CASE (T_deriv_test)
      {
        jif3D::ThreeDMagneticModel MagTest;
        const size_t nmeas = 3;
        const size_t ncells = 5;
        MakeRandomModel(MagTest, ncells, nmeas, false);
        double inclination = drand48();
        double declination = drand48();
        double fieldstrength = 1.0 + drand48();

        typedef typename jif3D::MinMemGravMagCalculator<jif3D::ThreeDMagneticModel> CalculatorType;
        boost::shared_ptr<jif3D::ThreeDGravMagImplementation<jif3D::ThreeDMagneticModel> > Implementation(
            new jif3D::OMPMagneticImp(inclination, declination, fieldstrength));

        boost::shared_ptr<CalculatorType> Calculator(new CalculatorType(Implementation));

        Calculator->SetDataTransform(
           boost::shared_ptr<jif3D::TotalField>(new jif3D::TotalField));

        jif3D::rvec Observed(Calculator->Calculate(MagTest));

        Observed *= 1.1;

        jif3D::ThreeDModelObjective<CalculatorType> Objective(*Calculator.get());
        Objective.SetObservedData(Observed);
        Objective.SetCoarseModelGeometry(MagTest);
        jif3D::rvec InvModel(MagTest.GetSusceptibilities().num_elements());
        std::copy(MagTest.GetSusceptibilities().origin(),
            MagTest.GetSusceptibilities().origin()
                + MagTest.GetSusceptibilities().num_elements(), InvModel.begin());

        CheckGradient(Objective, InvModel);
      }
    BOOST_AUTO_TEST_SUITE_END()
