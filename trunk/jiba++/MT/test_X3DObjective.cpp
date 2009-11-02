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
#include <boost/filesystem.hpp>
#include "X3DObjective.h"
#include "X3DModel.h"
#include "X3DMTCalculator.h"
#include "ReadWriteX3D.h"
#include "MTEquations.h"
#include "ReadWriteImpedances.h"

BOOST_AUTO_TEST_SUITE( X3DCalculator_Suite )

bool Between(const double limit1, const double limit2, const double value)
  {
    const double upper = std::max(limit1,limit2);
    const double lower = std::min(limit1,limit2);
    return (lower <= value) && (upper >= value);
  }

BOOST_AUTO_TEST_CASE  (X3D_basic_deriv_test)
    {
      const size_t xsize = 3;
      const size_t ysize = 6;
      const size_t zsize = 4;
      const size_t nbglayers = 5;
      const size_t nmod = xsize * ysize * zsize;
      jiba::X3DModel Model;

      Model.SetZCellSizes().resize(boost::extents[zsize]);

      Model.SetConductivities().resize(boost::extents[xsize][ysize][zsize]);
      std::vector<double> bg_thicknesses(nbglayers),bg_conductivities(nbglayers);

      const double deltax = 100.0;
      const double deltay = 100.0;
      const double deltaz = 100.0;
      Model.SetHorizontalCellSize(deltax,deltay,xsize,ysize);

      std::fill_n(Model.SetZCellSizes().origin(),zsize,deltaz);
      std::fill_n(Model.SetConductivities().origin(),nmod,0.01);
      std::fill_n(bg_conductivities.begin(),nbglayers,0.011);
      std::fill_n(bg_thicknesses.begin(),nbglayers,100.0);

      Model.SetBackgroundConductivities(bg_conductivities);
      Model.SetBackgroundThicknesses(bg_thicknesses);
      Model.SetFrequencies().push_back(1.0);
      Model.SetFrequencies().push_back(10.0);

      for (size_t i = 0; i < xsize; ++i)
        {
          for (size_t j = 0; j < ysize; ++j)
            {
              double currx = Model.GetXCoordinates()[i]+deltax/2.0;
              double curry = Model.GetYCoordinates()[j]+deltay/2.0;
              Model.AddMeasurementPoint(currx,curry,0.0);
            }
        }
      jiba::X3DModel TrueModel(Model);
      std::fill_n(TrueModel.SetConductivities().origin(),nmod,0.012);


      jiba::X3DMTCalculator Calculator;
      jiba::rvec Observed = Calculator.Calculate(TrueModel);


      std::vector<double> Freq(TrueModel.GetFrequencies());

      jiba::WriteImpedancesToNetCDF("gradimp.nc",Freq,TrueModel.GetMeasPosX(),TrueModel.GetMeasPosY(),TrueModel.GetMeasPosZ(),Observed);

      jiba::X3DObjective Objective;
      Objective.SetObservedData(Observed);
      Objective.SetModelGeometry(Model);
      jiba::rvec ModelVec(nmod);
      std::copy(Model.GetConductivities().origin(),Model.GetConductivities().origin()+nmod,ModelVec.begin());
      double misfit = Objective.CalcMisfit(ModelVec);
      BOOST_CHECK(misfit > 0.0);
      jiba::rvec Gradient = Objective.CalcGradient(ModelVec);


      /*std::ofstream outfile("grad.comp");

      for (size_t index = 0; index < nmod; ++index)
        {
          double delta = ModelVec(index) * 0.001;
          jiba::rvec Forward(ModelVec);
          jiba::rvec Backward(ModelVec);
          Forward(index) += delta;
          Backward(index) -= delta;
          double ForFDGrad = (Objective.CalcMisfit(Forward) - misfit)/(delta);
          double BackFDGrad = (misfit - Objective.CalcMisfit(Backward))/delta;
          BOOST_CHECK(Between(ForFDGrad,BackFDGrad,Gradient(index)));
          outfile << index << " " << ForFDGrad << " "<< BackFDGrad << " " << Gradient(index) << std::endl;
        }*/
    }

  BOOST_AUTO_TEST_SUITE_END()
