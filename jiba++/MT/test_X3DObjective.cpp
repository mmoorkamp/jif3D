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

#include "X3DObjective.h"
#include "X3DModel.h"
#include "X3DMTCalculator.h"
#include "MTEquations.h"

BOOST_AUTO_TEST_SUITE( X3DCalculator_Suite )

BOOST_AUTO_TEST_CASE  (X3D_basic_deriv_test)
    {
      const size_t xsize = 5;
      const size_t ysize = 5;
      const size_t zsize = 5;
      const size_t nbglayers = 5;
      const size_t nmod = xsize * ysize * zsize;
      jiba::X3DModel Model;

      Model.SetZCellSizes().resize(boost::extents[zsize]);

      Model.SetConductivities().resize(boost::extents[xsize][ysize][zsize]);
      std::vector<double> bg_thicknesses(nbglayers),bg_conductivities(nbglayers);

      const double deltax = 100.0;
      const double deltay = 100.0;
      const double deltaz = 100.0;
      Model.SetHorizontalCellSize(deltax,xsize,ysize);

      std::fill_n(Model.SetZCellSizes().origin(),zsize,deltaz);
      std::fill_n(Model.SetConductivities().origin(),nmod,0.01);
      std::fill_n(bg_conductivities.begin(),nbglayers,0.011);
      std::fill_n(bg_thicknesses.begin(),nbglayers,100.0);

      Model.SetBackgroundConductivities(bg_conductivities);
      Model.SetBackgroundThicknesses(bg_thicknesses);
      Model.SetFrequencies().push_back(10.0);

      jiba::X3DMTCalculator Calculator;
      jiba::rvec Impedance = Calculator.Calculate(Model);

      Impedance *= 1.01;
      jiba::X3DObjective Objective;
      Objective.SetObservedData(Impedance);
      Objective.SetModelGeometry(Model);
      jiba::rvec ModelVec(nmod);
      std::copy(Model.GetConductivities().origin(),Model.GetConductivities().origin()+nmod,ModelVec.begin());
      double misfit = Objective.CalcMisfit(ModelVec);
      BOOST_CHECK(misfit > 0.0);
      jiba::rvec Gradient = Objective.CalcGradient(ModelVec);

      //const size_t index =  62; //rand() % nmod;
      /*jiba::rvec FDGrad(nmod);
      std::ofstream outfile("grad.comp");
      for (size_t index = 0; index < nmod; ++index)
        {
          double delta = ModelVec(index) * 0.01;
          jiba::rvec Forward(ModelVec);
          jiba::rvec Backward(ModelVec);
          Forward(index) += delta;
          Backward(index) -= delta;
          FDGrad(index) = (Objective.CalcMisfit(Forward) - Objective.CalcMisfit(Backward))/(2*delta);
          outfile << index << " " << FDGrad(index) << " " << Gradient(index) << std::endl;
        }*/
      //std::cout << "Index: " << index << std::endl;
      //BOOST_CHECK_CLOSE(FDGrad,Gradient(index),1.0);
      // std::cout << "Gradient: " << Gradient << std::endl;
    }

  BOOST_AUTO_TEST_SUITE_END()
