//============================================================================
// Name        : test_inversion.h
// Author      : Max Moorkamp
// Version     :
// Copyright   : 2008, MM
//============================================================================

#define BOOST_TEST_MODULE Distributor test
#define BOOST_TEST_MAIN ...
#include <boost/test/included/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <stdlib.h>
#include "ModelDistributor.h"
#include "ModelTransforms.h"


//a general function to test various ModelTransforms
template<class TransformType>
void TestTransform()
  {
    //this number is arbitrary
    const size_t nelements = 11;
    jiba::rvec Physical(nelements), Reference(nelements),
        Generalized(nelements), Derivative(nelements), Ones(nelements);
    //some transforms assume the normal ranges for seismic velocities or densities
    //so we make the physical parameters have similar ranges
    //otherwise everything is random
    for (size_t i = 0; i < nelements; ++i)
      {
        Physical( i) = 2.0 + drand48();
        Reference( i) = 3.1 + drand48();
        Derivative( i) = drand48();
        Ones( i) = 1.0;
      }
    TransformType Transform(Reference);
    Generalized = Transform.PhysicalToGeneralized(Physical);
    jiba::rvec Compare = Transform.GeneralizedToPhysical(Generalized);
    //first we check whether backwards and forwards transformation produces
    //the right result
    for (size_t i = 0; i < nelements; ++i)
      {
        BOOST_CHECK_CLOSE(Physical(i),Compare(i),1e-5);
      }
    //then we check the derivative from the transformation class
    //against a finite difference
    for (size_t i = 0; i < nelements; ++i)
      {
        jiba::rvec PlusVec(Generalized), MinusVec(Generalized);
        const double delta = 0.001;
        const double h = Generalized(i) * delta;
        PlusVec(i) += h;
        MinusVec(i) -= h;
        const double DiffDeriv = (Transform.GeneralizedToPhysical(PlusVec)(i)
            - Transform.GeneralizedToPhysical(MinusVec)(i)) / (2 * h);
        jiba::rvec TransDeriv = Transform.Derivative(Generalized, Ones);
        BOOST_CHECK_CLOSE(DiffDeriv,TransDeriv(i),1e-2);
      }
    //finally we check whether the call to derivative applies
    //the chain rule properly
    jiba::rvec TransDeriv = Transform.Derivative(Generalized, Derivative);
    jiba::rvec CompDeriv = ublas::element_prod(Transform.Derivative(
        Generalized, Ones), Derivative);
    for (size_t i = 0; i < nelements; ++i)
      {
        BOOST_CHECK_CLOSE(TransDeriv(i),CompDeriv(i),1e-5);
      }
  }

BOOST_AUTO_TEST_SUITE( Distributor_Test_Suite )

BOOST_AUTO_TEST_CASE (basic_copy_test)
    {
      jiba::ModelDistributor Distributor;
      Distributor.AddTransformer(
          boost::shared_ptr<jiba::GeneralModelTransform>(
              new jiba::ModelCopyTransform()));
      Distributor.AddTransformer(
          boost::shared_ptr<jiba::GeneralModelTransform>(
              new jiba::ModelCopyTransform()));
      jiba::rvec Input(50);
      std::generate(Input.begin(), Input.end(), rand);

      jiba::rvec Output1(Distributor(Input, 0));
      jiba::rvec Output2(Distributor(Input, 1));

      BOOST_CHECK(std::equal(Input.begin(),Input.end(),Output1.begin()));
      BOOST_CHECK(std::equal(Input.begin(),Input.end(),Output2.begin()));
    }

  BOOST_AUTO_TEST_CASE (transform_test)
    {
      TestTransform<jiba::NormalizeTransform> ();
      TestTransform<jiba::LogTransform> ();
      TestTransform<jiba::LogDensityTransform> ();
      TestTransform<jiba::VelTransform> ();
      TestTransform<jiba::VelDensTransform> ();
    }
BOOST_AUTO_TEST_SUITE_END()
