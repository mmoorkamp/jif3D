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
#include <boost/test/test_tools.hpp>
#include <boost/filesystem.hpp>
#include "../Inversion/ThreeDModelObjective.h"
#include "X3DModel.h"
#include "X3DMTCalculator.h"
#include "ReadWriteX3D.h"
#include "MTEquations.h"
#include "ReadWriteImpedances.h"
#define BOOST_UBLAS_TYPE_CHECK_MIN (real_type(1))
BOOST_AUTO_TEST_SUITE( X3DObjective_Suite )

void  MakeMTModel(jif3D::X3DModel &Model)
    {
      const size_t xsize = 3;
      const size_t ysize = 4;
      const size_t zsize = 2;
      const size_t nbglayers = 5;
      const size_t nmod = xsize * ysize * zsize;

      Model.SetZCellSizes().resize(boost::extents[zsize]);

      Model.SetConductivities().resize(boost::extents[xsize][ysize][zsize]);
      std::vector<double> bg_thicknesses(nbglayers), bg_conductivities(nbglayers);

      const double deltax = 100.0;
      const double deltay = 100.0;
      const double deltaz = 100.0;
      Model.SetHorizontalCellSize(deltax, deltay, xsize, ysize);

      std::fill_n(Model.SetZCellSizes().origin(), zsize, deltaz);
      std::fill_n(Model.SetConductivities().origin(), nmod, 0.01);
      std::fill_n(bg_conductivities.begin(), nbglayers, 0.011);
      std::fill_n(bg_thicknesses.begin(), nbglayers, 100.0);

      Model.SetBackgroundConductivities(bg_conductivities);
      Model.SetBackgroundThicknesses(bg_thicknesses);
      Model.SetFrequencies().push_back(1.0);
      Model.SetFrequencies().push_back(2.0);
      Model.SetFrequencies().push_back(5.0);
      Model.SetFrequencies().push_back(10.0);
      srand48(time(0));
      for (size_t i = 1; i < xsize-1; ++i)
        {
          for (size_t j = 1; j < ysize-1; ++j)
            {
              double currx = Model.GetXCoordinates()[i] + deltax / 2.0;
              double curry = Model.GetYCoordinates()[j] + deltay / 2.0;
              double currz =  drand48() * deltaz * zsize;
              Model.AddMeasurementPoint(currx, curry, currz);
            }
        }
    }

  bool Between(const double limit1, const double limit2, const double value)
    {
      const double upper = std::max(limit1, limit2);
      const double lower = std::min(limit1, limit2);
      return (lower <= value) && (upper >= value);
    }

  BOOST_AUTO_TEST_CASE (X3D_fail_test)
    {
      jif3D::X3DModel Model;
      jif3D::rvec Observed;
      jif3D::X3DMTCalculator Calculator;
      jif3D::ThreeDModelObjective<jif3D::X3DMTCalculator> Objective(Calculator);
      //BOOST_CHECK_THROW(Objective.SetObservedData(Observed),jif3D::FatalException);
      //BOOST_CHECK_THROW(Objective.SetCoarseModelGeometry(Model),jif3D::FatalException);
      Observed.resize(10);
      Observed.clear();
      //BOOST_CHECK_NO_THROW(Objective.SetObservedData(Observed));
      MakeMTModel(Model);
      //BOOST_CHECK_NO_THROW(Objective.SetCoarseModelGeometry(Model));
      Model.ClearMeasurementPoints();
      Model.AddMeasurementPoint(10.0,12.0,0.0);
      Model.AddMeasurementPoint(13.0,14.0,30.0);
      Objective.SetCoarseModelGeometry(Model);
      //BOOST_CHECK_THROW(Objective.CalcMisfit(jif3D::rvec(Model.GetConductivities().num_elements())),jif3D::FatalException);
    }

  BOOST_AUTO_TEST_CASE (X3D_basic_deriv_test)
    {

      jif3D::X3DModel Model;
      MakeMTModel(Model);
      const size_t xsize = Model.GetXCoordinates().size();
      const size_t ysize = Model.GetYCoordinates().size();
      const size_t zsize = Model.GetZCoordinates().size();
      const size_t nmod = xsize * ysize * zsize;

      jif3D::X3DModel TrueModel(Model);
      std::fill_n(TrueModel.SetConductivities().origin(),nmod,0.012);

      jif3D::X3DMTCalculator Calculator;
      jif3D::rvec Observed = Calculator.Calculate(TrueModel);

      std::vector<double> Freq(TrueModel.GetFrequencies());

      jif3D::WriteImpedancesToNetCDF("gradimp.nc",Freq,TrueModel.GetMeasPosX(),TrueModel.GetMeasPosY(),TrueModel.GetMeasPosZ(),Observed);

      jif3D::ThreeDModelObjective<jif3D::X3DMTCalculator> Objective(Calculator);
      Objective.SetObservedData(Observed);
      Objective.SetCoarseModelGeometry(Model);
      jif3D::rvec ModelVec(nmod);
      std::copy(Model.GetConductivities().origin(),Model.GetConductivities().origin()+nmod,ModelVec.begin());
      double misfit = Objective.CalcMisfit(ModelVec);
      BOOST_CHECK(misfit > 0.0);
      jif3D::rvec Gradient = Objective.CalcGradient(ModelVec);

      std::ofstream outfile("grad.comp");

      for (size_t index = 0; index < nmod; ++index)
        {
          double delta = ModelVec(index) * 0.001;
          jif3D::rvec Forward(ModelVec);
          jif3D::rvec Backward(ModelVec);
          Forward(index) += delta;
          Backward(index) -= delta;
          double ForFDGrad = (Objective.CalcMisfit(Forward) - misfit)/(delta);
          double BackFDGrad = (misfit - Objective.CalcMisfit(Backward))/delta;
          double CentFDGrad = (ForFDGrad + BackFDGrad)/2.0;
          BOOST_CHECK(Between(ForFDGrad,BackFDGrad,Gradient(index)) || (CentFDGrad - Gradient(index))/CentFDGrad < 0.01);
          outfile << index << " " << ForFDGrad << " "<< BackFDGrad << " " << (ForFDGrad + BackFDGrad)/2.0 << " " << Gradient(index) << std::endl;
        }
    }


  BOOST_AUTO_TEST_SUITE_END()
