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
#include "X3DObjective.h"
#include "X3DModel.h"
#include "X3DMTCalculator.h"
#include "ReadWriteX3D.h"
#include "MTEquations.h"
#include "ReadWriteImpedances.h"

BOOST_AUTO_TEST_SUITE( X3DCalculator_Suite )

void  MakeMTModel(jiba::X3DModel &Model)
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

      for (size_t i = 0; i < xsize; ++i)
        {
          for (size_t j = 0; j < ysize; ++j)
            {
              double currx = Model.GetXCoordinates()[i] + deltax / 2.0;
              double curry = Model.GetYCoordinates()[j] + deltay / 2.0;
              Model.AddMeasurementPoint(currx, curry, 250.0);
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
      jiba::X3DModel Model;
      jiba::rvec Observed;
      jiba::X3DObjective Objective;
      BOOST_CHECK_THROW(Objective.SetObservedData(Observed),jiba::FatalException);
      BOOST_CHECK_THROW(Objective.SetModelGeometry(Model),jiba::FatalException);
      Observed.resize(10);
      Observed.clear();
      BOOST_CHECK_NO_THROW(Objective.SetObservedData(Observed));
      MakeMTModel(Model);
      BOOST_CHECK_NO_THROW(Objective.SetModelGeometry(Model));
      Model.ClearMeasurementPoints();
      Model.AddMeasurementPoint(10.0,12.0,250.0);
      Model.AddMeasurementPoint(13.0,14.0,350.0);
      Objective.SetModelGeometry(Model);
      BOOST_CHECK_THROW(Objective.CalcMisfit(jiba::rvec(Model.GetConductivities().num_elements())),jiba::FatalException);
    }

  BOOST_AUTO_TEST_CASE (X3D_basic_deriv_test)
    {

      jiba::X3DModel Model;
      MakeMTModel(Model);
      const size_t xsize = Model.GetXCoordinates().size();
      const size_t ysize = Model.GetYCoordinates().size();
      const size_t zsize = Model.GetZCoordinates().size();
      const size_t nmod = xsize * ysize * zsize;

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

      std::ofstream outfile("grad.comp");

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
        }
    }

  BOOST_AUTO_TEST_CASE (X3D_gradient_reproc_test)
    {
      jiba::X3DModel Model;
      MakeMTModel(Model);
      const size_t xsize = Model.GetXCoordinates().size();
      const size_t ysize = Model.GetYCoordinates().size();
      const size_t zsize = Model.GetZCoordinates().size();
      const size_t nmod = xsize * ysize * zsize;
      const size_t nruns = 5;

      jiba::X3DModel TrueModel(Model);
      std::fill_n(TrueModel.SetConductivities().origin(),nmod,0.1);

      jiba::X3DMTCalculator Calculator;
      jiba::rvec Observed = Calculator.Calculate(TrueModel);
      std::vector<double> FitVec;
      std::vector<jiba::rvec> GradVec;
      std::vector<jiba::rvec> DiffVec;
      std::vector<boost::shared_ptr<jiba::X3DObjective> > Objectives;
      for (size_t i = 0; i < nruns; ++i)
        {
          boost::shared_ptr<jiba::X3DObjective> Objective(new jiba::X3DObjective);
          Objectives.push_back(Objective);
          Objective->SetObservedData(Observed);
          Objective->SetModelGeometry(Model);
          jiba::rvec ModelVec(nmod);
          std::copy(Model.GetConductivities().origin(),Model.GetConductivities().origin()+nmod,ModelVec.begin());
          double misfit = Objective->CalcMisfit(ModelVec);
          DiffVec.push_back(Objective->GetDataDifference());
          jiba::rvec Grad = Objective->CalcGradient(ModelVec);
          FitVec.push_back(misfit);
          GradVec.push_back(Grad);
        }
      for (size_t i = 1; i < nruns; ++i)
        {
          BOOST_CHECK(FitVec.at(0) == FitVec.at(i));
          for (size_t j = 0; j < GradVec.size(); ++j)
            {
              BOOST_CHECK_CLOSE(GradVec.at(0)(j),GradVec.at(i)(j),std::numeric_limits<double>::epsilon());
              if (GradVec.at(0)(j) != GradVec.at(i)(j))
                {
                  std::cout << "i: " << i << " j: " << j << std::endl;
                }
            }

          for (size_t j = 0; j < DiffVec.size(); ++j)
            {
              BOOST_CHECK(DiffVec.at(0)(j) == DiffVec.at(i)(j));
            }

        }
    }
  BOOST_AUTO_TEST_SUITE_END()
