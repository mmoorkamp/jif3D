//============================================================================
// Name        : test_transform.cpp
// Author      : Max Moorkamp
// Version     :
// Copyright   : 2008, MM
//============================================================================

#define BOOST_TEST_MODULE Transform test
#define BOOST_TEST_MAIN ...
#include "../Global/Jif3DTesting.h"
#include <boost/test/included/unit_test.hpp>
#include <boost/test/tools/floating_point_comparison.hpp>
#include <boost/lambda/lambda.hpp>
#include <boost/math/special_functions/fpclassify.hpp>
#include <stdlib.h>
#include "../Inversion/ModelTransforms.h"
#include "../Gravity/ThreeDGravityModel.h"
#include "../Global/Jif3DPlatformHelper.h"

/*!
 * \file test_transform.cpp
 * This file contains some tests necessary to verify the functionality of the
 * different model transforms.
 */

//a general function to test various ModelTransforms
//it checks whether consecutive forward and inverse transforms
//give back the original vector and checks the implementation
//of the gradient against a finite difference approximation
void TestTransform(const jif3D::GeneralModelTransform &Transform, const size_t nelements,
    const double min = 2.0, const double max = 3.0)
  {
    jif3D::platform::srand48((unsigned int)time(NULL));
    jif3D::rvec Physical(nelements), Generalized(nelements), Derivative(nelements), One(
        nelements);
    //some transforms assume the normal ranges for seismic velocities or densities
    //so we make the physical parameters have similar ranges
    //otherwise everything is random
    for (size_t i = 0; i < nelements; ++i)
      {
        Physical(i) = min + (max - min) * jif3D::platform::drand48();
        Derivative(i) = drand48() - 1.0;
      }
    //first we check whether backwards and forwards transformation produces
    //the right result
    Generalized = Transform.PhysicalToGeneralized(Physical);
    jif3D::rvec Compare = Transform.GeneralizedToPhysical(Generalized);
    for (size_t i = 0; i < nelements; ++i)
      {
        BOOST_CHECK_CLOSE(Physical(i), Compare(i), 1e-5);
      }
    //then we check the derivative from the transformation class
    //against a finite difference
    for (size_t i = 0; i < nelements; ++i)
      {
        One.clear();
        One(i) = 1.0;
        jif3D::rvec TransDeriv = Transform.Derivative(Generalized, One);
        jif3D::rvec PlusVec(Generalized), MinusVec(Generalized);
        const double delta = 0.001;
        const double h =  delta;
        PlusVec(i) += h;
        MinusVec(i) -= h;
        const double DiffDeriv = (Transform.GeneralizedToPhysical(PlusVec)(i)
            - Transform.GeneralizedToPhysical(MinusVec)(i)) / (2 * h);
        BOOST_CHECK_CLOSE(DiffDeriv, TransDeriv(i), 0.5);
      }

  }

BOOST_AUTO_TEST_SUITE (Transform_Test_Suite)
//Test the functionality of ModelCopyTransform
//This transformation simply copies its input arguments
    BOOST_AUTO_TEST_CASE (basic_copy_test)
      {
        jif3D::ModelCopyTransform CTrans;
        //generate a random input vector
        const size_t nelements = 5 + rand() % 100;
        jif3D::rvec Input(nelements);
        std::generate(Input.begin(), Input.end(), jif3D::platform::drand48);
        //for this particular transform each transformation
        //direction should always return the original result
        jif3D::rvec Output1(CTrans.GeneralizedToPhysical(Input));
        jif3D::rvec Output2(CTrans.PhysicalToGeneralized(Input));
        //check that all values are identical
        BOOST_CHECK(std::equal(Input.begin(), Input.end(), Output1.begin()));
        BOOST_CHECK(std::equal(Input.begin(), Input.end(), Output2.begin()));
        TestTransform(CTrans, nelements);
      }

    BOOST_AUTO_TEST_CASE (Normalize_transform_test)
      {
        //generate a random vector
        const size_t nelements = 5 + rand() % 100;
        jif3D::rvec Reference(nelements);
        //we need the elements of the reference vector
        //to be > 0, so we have to add a number to the
        //output of the random generator, this is a pain
        //with generate and lambda function, so we use a plain
        //old loop
        for (size_t i = 0; i < nelements; ++i)
          {
            Reference(i) = 2.0 + jif3D::platform::drand48();
          }
        TestTransform(jif3D::NormalizeTransform(Reference), nelements);
      }

    BOOST_AUTO_TEST_CASE (Log_transform_test)
      {

        const size_t nelements = 5 + rand() % 100;
        jif3D::rvec Reference(nelements);
        for (size_t i = 0; i < nelements; ++i)
          {
            Reference(i) = 2.0 + jif3D::platform::drand48();
          }
        TestTransform(jif3D::LogTransform(Reference), nelements);
      }

    BOOST_AUTO_TEST_CASE (Wavelet_transform_test)
      {
        //for the wavelet transform we need
        //to specify the 3 dimensions of the grid
        //each dimension has to be a power of 2
        const size_t nx = 4;
        const size_t ny = 8;
        const size_t nz = 8;
        const size_t nelements = nx * ny * nz;

        jif3D::WaveletModelTransform WT(nx, ny, nz);
        TestTransform(WT, nelements);
      }

    BOOST_AUTO_TEST_CASE (wavelet_throw_test)
      {
        //check that WaveletModelTransform throws when one dimension is not a power of two
        BOOST_CHECK_THROW(jif3D::WaveletModelTransform(3, 4, 8), jif3D::FatalException);
        BOOST_CHECK_THROW(jif3D::WaveletModelTransform(4, 5, 8), jif3D::FatalException);
        BOOST_CHECK_THROW(jif3D::WaveletModelTransform(4, 4, 9), jif3D::FatalException);
      }

    BOOST_AUTO_TEST_CASE (Tanh_transform_test)
      {
        const size_t nelements = 5 + rand() % 100;
        TestTransform(jif3D::TanhTransform(0.0, 1000), nelements);
      }

    BOOST_AUTO_TEST_CASE (Tanh_vector_transform_test)
      {
        const size_t nelements = 5 + rand() % 100;
        jif3D::rvec min(nelements), max(nelements);
        std::generate(min.begin(),min.end(),[](){return drand48();});
        std::generate(max.begin(),max.end(),[](){return 3.0 + 5 * drand48();});
        TestTransform(jif3D::TanhTransform(min, max), nelements);
      }

    BOOST_AUTO_TEST_CASE (LogLim_transform_test)
      {
        const size_t nelements = 5 + rand() % 100;
        TestTransform(jif3D::LogLimTransform(0.01, 1000), nelements);
      }

    BOOST_AUTO_TEST_CASE (DensPVel_transform_test)
      {
        const size_t nelements = (5 + rand() % 100) * 3;
        std::vector<double> bg_dens(nelements/3);
        std::generate(bg_dens.begin(),bg_dens.end(),[](){return drand48();});
        //current testing is not made for this transform that does not preserve number of elements
       // TestTransform(jif3D::DensPVelTransform(bg_dens), nelements);
      }

    BOOST_AUTO_TEST_CASE (Density_transform_test)
      {
        //for density transform we can use a model object
        //to indicate where we want to apply the relationship
        //so we use a ThreeDGravityModel object
        const size_t nx = 5;
        const size_t ny = 7;
        const size_t nz = 9;
        const size_t nelements = nx * ny * nz;

        jif3D::ThreeDGravityModel Model;
        Model.SetDensities().resize(boost::extents[nx][ny][nz]);
        std::fill_n(Model.SetDensities().origin(), nelements, 1.0);
        //the density transform needs a transform that translates
        //the general parameters into slowness, before it then
        //translates slowness to density. Here we use TanhTransform
        //which does not make physical sense, but it is non-linear
        // and therefore well suited for testing the combined
        //transformation and gradient calculation
        TestTransform(
            jif3D::DensityTransform(
                boost::shared_ptr<jif3D::GeneralModelTransform>(
                    new jif3D::TanhTransform(0.0, 1000.0)), Model, 1.0), nelements);
      }

    BOOST_AUTO_TEST_CASE (Conductivity_transform_test)
      {
        //for conductivity transform we can use a model object
        //to indicate where we want to apply the relationship
        //so we use a ThreeDGravityModel object
        const size_t nx = 4;
        const size_t ny = 8;
        const size_t nz = 8;
        const size_t nelements = nx * ny * nz;

        jif3D::ThreeDGravityModel Model;
        Model.SetDensities().resize(boost::extents[nx][ny][nz]);
        std::fill_n(Model.SetDensities().origin(), nelements, 1.0);
        //the conductivity transform needs a transform that translates
        //the general parameters into slowness, before it then
        //translates slowness to conductivity. Here we use TanhTransform
        //which does not make physical sense, but it is non-linear
        // and therefore well suited for testing the combined
        //transformation and gradient calculation
        TestTransform(
            jif3D::ConductivityTransform(
                boost::shared_ptr<jif3D::GeneralModelTransform>(
                    new jif3D::TanhTransform(0.0, 10.0)), Model, 1.0), nelements, 0.001,
            0.5);
      }

    BOOST_AUTO_TEST_CASE (FChained_transform_test)
      {
        //test the chaining of several transformations
        //here we append consecutive transforms
        //in the order they will be applied
        const size_t nelements = 5 + rand() % 100;
        jif3D::rvec Reference(nelements);
        for (size_t i = 0; i < nelements; ++i)
          {
            Reference(i) = 2.0 + jif3D::platform::drand48();
          }

        jif3D::ChainedTransform TransformForward;
        TransformForward.AppendTransform(
            boost::shared_ptr<jif3D::GeneralModelTransform>(
                new jif3D::TanhTransform(-5, 5.0)));
        TransformForward.AppendTransform(
            boost::shared_ptr<jif3D::GeneralModelTransform>(
                new jif3D::LogTransform(Reference)));
        TestTransform(TransformForward, nelements);
      }

    BOOST_AUTO_TEST_CASE (RCHained_transform_test)
      {
        const size_t nelements = 5 + rand() % 100;
        jif3D::rvec Reference(nelements);
        for (size_t i = 0; i < nelements; ++i)
          {
            Reference(i) = 2.0 + jif3D::platform::drand48();
          }
        //test the chaining of several transformations
        //here we prepend consecutive transforms, i.e.
        //we specify them in the reverse order they will be applied
        jif3D::ChainedTransform TransformReverse;
        TransformReverse.PrependTransform(
            boost::shared_ptr<jif3D::GeneralModelTransform>(
                new jif3D::LogTransform(Reference)));
        TransformReverse.PrependTransform(
            boost::shared_ptr<jif3D::GeneralModelTransform>(
                new jif3D::TanhTransform(-5, 5.0)));
        TestTransform(TransformReverse, nelements);
      }
    BOOST_AUTO_TEST_SUITE_END()
