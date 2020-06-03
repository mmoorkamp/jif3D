//============================================================================
// Name        : test_X3DModel.cpp
// Author      : 8 Mar 2018
// Version     : 
// Copyright   : 2018, mm489
//============================================================================

#define BOOST_TEST_MODULE X3DModel test
#define BOOST_TEST_MAIN ...
#include "../Global/Jif3DTesting.h"
#include <boost/test/included/unit_test.hpp>
#include <boost/test/tools/floating_point_comparison.hpp>
#include "X3DModel.h"
#include "../Global/Jif3DPlatformHelper.h"

BOOST_AUTO_TEST_SUITE( X3DModel_Suite )

    BOOST_AUTO_TEST_CASE (X3D_Comparison_test)
      {
        jif3D::X3DModel TestModel;
        const size_t xsize = rand() % 10 + 2;
        const size_t ysize = rand() % 10 + 2;
        const size_t zsize = rand() % 10 + 2;
        const size_t nbglayers = rand() % 10 + 2;
        const double deltax = 11.2;
        const double deltay = 17.4;
        TestModel.SetMeshSize(xsize, ysize, zsize);
        TestModel.SetHorizontalCellSize(deltax, deltay, xsize, ysize);
        std::vector<double> Depths(zsize);
        std::generate(Depths.begin(), Depths.end(), jif3D::platform::drand48);
        TestModel.SetZCellSizes() = Depths;

        jif3D::ThreeDModelBase::t3DModelData Data(boost::extents[xsize][ysize][zsize]);
        std::vector<double> bg_thicknesses(nbglayers), bg_conductivities(nbglayers);

        std::generate_n(TestModel.SetZCellSizes().begin(), zsize,
            jif3D::platform::drand48);
        std::generate_n(TestModel.SetData().origin(), xsize * ysize * zsize,
            jif3D::platform::drand48);
        std::generate_n(bg_thicknesses.begin(), nbglayers, jif3D::platform::drand48);
        std::generate_n(bg_conductivities.begin(), nbglayers, jif3D::platform::drand48);
        TestModel.SetBackgroundThicknesses(bg_thicknesses);
        TestModel.SetBackgroundConductivities(bg_conductivities);
        jif3D::X3DModel CompareModel(TestModel);
        BOOST_CHECK(TestModel == CompareModel);
        CompareModel.SetMeshSize(xsize + 1, ysize, zsize);
        BOOST_CHECK(!(TestModel == CompareModel));

        CompareModel = TestModel;
        BOOST_CHECK(TestModel == CompareModel);
        CompareModel.SetData()[0][0][0] *= 1.1;
        BOOST_CHECK(!(TestModel == CompareModel));

        CompareModel = TestModel;
        CompareModel.AddMeasurementPoint(0, 1, 2);
        BOOST_CHECK(!(TestModel == CompareModel));

        CompareModel = TestModel;
        CompareModel.SetOrigin(1, 2, 3);
        BOOST_CHECK(!(TestModel == CompareModel));

        CompareModel = TestModel;
        CompareModel.SetHorizontalCellSize(deltax + 1, deltay, xsize, ysize);
        BOOST_CHECK(!(TestModel == CompareModel));

        CompareModel = TestModel;
        CompareModel.SetHorizontalCellSize(deltax, deltay + 1.0, xsize, ysize);
        BOOST_CHECK(!(TestModel == CompareModel));

        CompareModel = TestModel;
        CompareModel.SetZCellSizes()[0] *= 1.1;
        BOOST_CHECK(!(TestModel == CompareModel));

        CompareModel = TestModel;
        bg_thicknesses[0] *= 1.1;
        CompareModel.SetBackgroundThicknesses(bg_thicknesses);
        BOOST_CHECK(!(TestModel == CompareModel));

        CompareModel = TestModel;
        bg_conductivities[0] *= 1.1;
        CompareModel.SetBackgroundConductivities(bg_conductivities);
        BOOST_CHECK(!(TestModel == CompareModel));
      }

    BOOST_AUTO_TEST_SUITE_END()
