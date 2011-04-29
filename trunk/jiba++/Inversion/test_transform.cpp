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
#include "WaveletModelTransform.h"

//a general function to test various ModelTransforms
void TestTransform(const jiba::GeneralModelTransform &Transform,
    const size_t nelements, const double min = 2.0, const double max = 3.0)
  {
    srand48(time(NULL));
    jiba::rvec Physical(nelements), Generalized(nelements), Derivative(
        nelements), One(nelements);
    //some transforms assume the normal ranges for seismic velocities or densities
    //so we make the physical parameters have similar ranges
    //otherwise everything is random
    for (size_t i = 0; i < nelements; ++i)
      {
        Physical(i) = min + (max - min) * drand48();
        Derivative(i) = drand48() - 1.0;
      }

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
        One.clear();
        One(i) = 1.0;
        jiba::rvec TransDeriv = Transform.Derivative(Generalized, One);
        jiba::rvec PlusVec(Generalized), MinusVec(Generalized);
        const double delta = 0.001;
        const double h = Generalized(i) * delta;
        PlusVec(i) += h;
        MinusVec(i) -= h;
        const double DiffDeriv = (Transform.GeneralizedToPhysical(PlusVec)(i)
            - Transform.GeneralizedToPhysical(MinusVec)(i)) / (2 * h);
        BOOST_CHECK_CLOSE(DiffDeriv,TransDeriv(i),0.1);
        //finally we check whether the call to derivative applies
        //the chain rule properly
        jiba::rvec ChainDeriv = Transform.Derivative(Generalized, Derivative);
        jiba::rvec CompDeriv = ublas::element_prod(Transform.Derivative(
            Generalized, One), Derivative);
        BOOST_CHECK_CLOSE(ChainDeriv(i),CompDeriv(i),1e-5);
      }

  }

BOOST_AUTO_TEST_SUITE( Distributor_Test_Suite )

BOOST_AUTO_TEST_CASE  (basic_copy_test)
    {
      srand48(time(NULL));
      jiba::ModelDistributor Distributor;
      Distributor.AddTransformer(
          boost::shared_ptr<jiba::GeneralModelTransform>(
              new jiba::ModelCopyTransform()));
      Distributor.AddTransformer(
          boost::shared_ptr<jiba::GeneralModelTransform>(
              new jiba::ModelCopyTransform()));
      jiba::rvec Input(50);
      std::generate(Input.begin(), Input.end(), drand48);

      jiba::rvec Output1(Distributor(Input, 0));
      jiba::rvec Output2(Distributor(Input, 1));

      BOOST_CHECK(std::equal(Input.begin(),Input.end(),Output1.begin()));
      BOOST_CHECK(std::equal(Input.begin(),Input.end(),Output2.begin()));
    }

  BOOST_AUTO_TEST_CASE (transform_test)
    {
      //for WaveletModelTransform we need all dimensions to be a power of two
      const size_t nx = 4;
      const size_t ny = 8;
      const size_t nz = 8;
      const size_t nelements = nx * ny *nz;
      jiba::rvec Reference(nelements);
      for (size_t i = 0; i < nelements; ++i)
        {
          Reference( i) = 2.0 + drand48();
        }
      TestTransform(jiba::NormalizeTransform(Reference),nelements);
      TestTransform(jiba::LogTransform(Reference),nelements);
      TestTransform(jiba::DensityTransform(boost::shared_ptr<jiba::GeneralModelTransform>(new jiba::LogTransform(Reference))),nelements);
      jiba::WaveletModelTransform WT(nx,ny,nz);
      std::cout << "Testing wavelet \n\n" << std::endl;
      TestTransform(WT,nelements);
      std::cout << "End Testing wavelet \n\n" << std::endl;
      TestTransform(jiba::TanhTransform(0.0,1000),nelements);
      TestTransform(jiba::DensityTransform(boost::shared_ptr<jiba::GeneralModelTransform>(new jiba::TanhTransform(0.0,1000.0))),nelements);
      TestTransform(jiba::ConductivityTransform(boost::shared_ptr<jiba::GeneralModelTransform>(new jiba::TanhTransform(0.0,10.0))),nelements,0.001,0.5);
      jiba::ChainedTransform Transform;
      Transform.AddTransform(boost::shared_ptr<jiba::GeneralModelTransform>(new jiba::TanhTransform(-5,5.0)) );
      Transform.AddTransform(boost::shared_ptr<jiba::GeneralModelTransform>(new jiba::LogTransform(Reference)));

      TestTransform(Transform,nelements);
    }

  BOOST_AUTO_TEST_CASE (wavelet_check_test)
    {
      //check that WaveletModelTransform throws when one dimension is not a power of two
      BOOST_CHECK_THROW(jiba::WaveletModelTransform(3,4,8),jiba::FatalException);
      BOOST_CHECK_THROW(jiba::WaveletModelTransform(4,5,8),jiba::FatalException);
      BOOST_CHECK_THROW(jiba::WaveletModelTransform(4,4,9),jiba::FatalException);
    }
  BOOST_AUTO_TEST_SUITE_END()
