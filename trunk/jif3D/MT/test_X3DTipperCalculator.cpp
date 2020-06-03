//============================================================================
// Name        : test_ReadWriteX3D.cpp
// Author      : Feb 17, 2009
// Version     :
// Copyright   : 2009, mmoorkamp
//============================================================================

#define BOOST_TEST_MODULE X3DCalculator test
#define BOOST_TEST_MAIN ...
#include "../Global/Jif3DTesting.h"
#include <boost/test/included/unit_test.hpp>
#include <boost/test/tools/floating_point_comparison.hpp>

#include "../Global/NumUtil.h"
#include "X3DTipperCalculator.h"
#include "X3DModel.h"
#include "ReadWriteImpedances.h"
#include "TipperData.h"

BOOST_AUTO_TEST_SUITE (X3DCalculator_Suite)

BOOST_AUTO_TEST_CASE (X3D_forward_hs_test)
  {
    //create a 3D version of a 1D model
    const size_t xsize = 20;
    const size_t ysize = 20;
    const size_t zsize = 10;
    const size_t nbglayers = 5;
    jif3D::X3DModel Model;
    jif3D::TipperData Data;
    Model.SetConductivities().resize(boost::extents[xsize][ysize][zsize]);
    std::vector<double> bg_thicknesses(nbglayers), bg_conductivities(nbglayers);

    const double deltax = 100.0;
    const double deltay = 100.0;
    const double deltaz = 100.0;
    const double freq = 1.0;
    const double cond = 0.01;
    Model.SetHorizontalCellSize(deltax, deltay, xsize, ysize);
    jif3D::ThreeDModelBase::t3DModelDim ZCS(zsize, deltaz);
    Model.SetZCellSizes(ZCS);
    std::fill_n(Model.SetConductivities().origin(), xsize * ysize * zsize, cond);
    std::fill_n(bg_conductivities.begin(), nbglayers, cond);
    bg_conductivities.back() *= 1.001;
    std::fill_n(bg_thicknesses.begin(), nbglayers, 200.0);

    Model.SetBackgroundConductivities(bg_conductivities);
    Model.SetBackgroundThicknesses(bg_thicknesses);
    Data.SetFrequencies(
          { freq});

    jif3D::X3DTipperCalculator Calculator;
    for (size_t i = 1; i < xsize / 2; ++i)
    for (size_t j = 1; j < ysize / 2; ++j)
      {
        Data.AddMeasurementPoint(Model.GetXCoordinates()[i] + deltax / 2.0,
            Model.GetYCoordinates()[j] + deltay / 2.0, 0.0);
      }
    std::vector<double> Dummy(Data.GetMeasPosX().size()*4);
    Data.SetDataAndErrors(Dummy,Dummy);
    Data.CompleteObject();
    jif3D::rvec Tipper = Calculator.Calculate(Model, Data);
    const double prec = 1e-5;
    for (size_t i = 0; i < Tipper.size() / 4; ++i)
      {
        BOOST_CHECK_SMALL(Tipper(i * 4 ), prec);
        BOOST_CHECK_SMALL(Tipper(i * 4 + 1), prec);
        BOOST_CHECK_SMALL(Tipper(i * 4 + 2), prec);
        BOOST_CHECK_SMALL(Tipper(i * 4 + 3), prec);
      }

  }

BOOST_AUTO_TEST_SUITE_END()
