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
#include <boost/test/floating_point_comparison.hpp>

#include "../Global/NumUtil.h"
#include "X3DTipperCalculator.h"
#include "X3DModel.h"
#include "ReadWriteImpedances.h"

BOOST_AUTO_TEST_SUITE( X3DCalculator_Suite )

    BOOST_AUTO_TEST_CASE (X3D_forward_hs_test)
      {
        //create a 3D version of a 1D model
        const size_t xsize = 20;
        const size_t ysize = 20;
        const size_t zsize = 10;
        const size_t nbglayers = 5;
        jif3D::X3DModel Model;

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
        Model.SetFrequencies().push_back(freq);

        jif3D::X3DTipperCalculator Calculator;
        for (size_t i = 1; i < xsize / 2; ++i)
          for (size_t j = 1; j < ysize / 2; ++j)
            {
              Model.AddMeasurementPoint(Model.GetXCoordinates()[i] + deltax / 2.0,
                  Model.GetYCoordinates()[j] + deltay / 2.0, 0.0);
            }

        jif3D::rvec Impedance = Calculator.Calculate(Model);
        const double prec = 0.05;
        for (size_t i = 0; i < Impedance.size() / 4; ++i)
          {
            BOOST_CHECK_CLOSE(Impedance(i * 4 ), 0.0, prec);
            BOOST_CHECK_CLOSE(Impedance(i * 4 + 1), 0.0, prec);
            BOOST_CHECK_CLOSE(Impedance(i * 4 + 2), 0.0, prec);
            BOOST_CHECK_CLOSE(Impedance(i * 4 + 3), 0.0, prec);
          }

      }


    BOOST_AUTO_TEST_SUITE_END()
