//============================================================================
// Name        : test_ReadWriteX3D.cpp
// Author      : Feb 17, 2009
// Version     :
// Copyright   : 2009, mmoorkamp
//============================================================================


#define BOOST_TEST_MODULE X3DCalculator test
#define BOOST_TEST_MAIN ...
#include <boost/test/included/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

#include "X3DMTCalculator.h"
#include "X3DModel.h"
#include "MTEquations.h"

BOOST_AUTO_TEST_SUITE( X3DCalculator_Suite )

BOOST_AUTO_TEST_CASE  (X3D_forward_hs_test)
    {
      //create a random number of cells and background layers
      const size_t xsize = 10;
      const size_t ysize = 10;
      const size_t zsize = 10;
      const size_t nbglayers = 5;
      jiba::X3DModel Model;

      Model.SetZCellSizes().resize(boost::extents[zsize]);

      Model.SetConductivities().resize(boost::extents[xsize][ysize][zsize]);
      std::vector<double> bg_thicknesses(nbglayers),bg_conductivities(nbglayers);

      const double deltax = 100.0;
      const double deltay = 100.0;
      const double deltaz = 100.0;
      Model.SetHorizontalCellSize(deltax,xsize,ysize);

      std::fill_n(Model.SetZCellSizes().origin(),zsize,deltaz);
      std::fill_n(Model.SetConductivities().origin(),xsize*ysize*zsize,100.0);
      std::fill_n(bg_conductivities.begin(),nbglayers,100.0);
      std::fill_n(bg_thicknesses.begin(),nbglayers,100.0);

      Model.SetBackgroundConductivities(bg_conductivities);
      Model.SetBackgroundThicknesses(bg_thicknesses);
      Model.SetFrequencies().push_back(10.0);

      jiba::X3DMTCalculator Calculator;
      jiba::rvec Impedance = Calculator.Calculate(Model);
      std::cout << Impedance << std::endl;
      std::complex<double> HsImp = jiba::ImpedanceHalfspace(10.0,100.0);
      BOOST_CHECK_CLOSE(Impedance(2),HsImp.real(),0.01);
      BOOST_CHECK_CLOSE(Impedance(3),HsImp.imag(),0.01);

    }

  BOOST_AUTO_TEST_SUITE_END()
