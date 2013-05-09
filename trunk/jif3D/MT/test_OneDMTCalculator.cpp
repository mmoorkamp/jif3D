//============================================================================
// Name        : test_ReadWriteX3D.cpp
// Author      : Feb 17, 2009
// Version     :
// Copyright   : 2009, mmoorkamp
//============================================================================

#define BOOST_TEST_MODULE OneDMTCalculator test
#define BOOST_TEST_MAIN ...
#include <boost/test/included/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

#include "../Global/Noise.h"
#include "../Global/NumUtil.h"
#include "OneDMTCalculator.h"
#include "OneDMTObjective.h"
#include "MTEquations.h"
#include "X3DModel.h"

BOOST_AUTO_TEST_SUITE( OneDMTCalculator_Suite )

bool  Between(const double limit1, const double limit2, const double value)
    {
      const double upper = std::max(limit1, limit2);
      const double lower = std::min(limit1, limit2);
      return (lower <= value) && (upper >= value);
    }

  BOOST_AUTO_TEST_CASE (OneD_forward_hs_test)
    {
      const size_t nbglayers = 5;
      const size_t nfreq = 6;
      jiba::X3DModel Model;
      std::vector<double> bg_thicknesses(nbglayers),bg_conductivities(nbglayers);

      const double freq = 1.0;
      const double cond = 0.01;

      std::fill_n(bg_conductivities.begin(),nbglayers,cond);
      std::fill_n(bg_thicknesses.begin(),nbglayers,200.0);

      Model.SetBackgroundConductivities(bg_conductivities);
      Model.SetBackgroundThicknesses(bg_thicknesses);
      for (size_t i = 0; i < nfreq; ++i)
        {
          Model.SetFrequencies().push_back(freq*(i+1));
        }

      jiba::OneDMTCalculator Calculator;

      jiba::rvec Impedance = Calculator.Calculate(Model);

      const double prec = 0.05;
      for (size_t i = 0; i < nfreq; ++i)
        {
          std::complex<double> HsImp = jiba::ImpedanceHalfspace(freq*(i+1),cond);
          BOOST_CHECK_CLOSE(Impedance(i*2),HsImp.real(),prec);
          BOOST_CHECK_CLOSE(Impedance(i*2+1),HsImp.imag(),prec);
        }

    }

  BOOST_AUTO_TEST_CASE (OneD_twolayer_test)
    {
      const size_t nbglayers = 5;

      jiba::X3DModel Model;
      std::vector<double> bg_thicknesses(nbglayers),bg_conductivities(nbglayers);

      const double freqhigh = 100.0;
      const double freqlow = 0.00001;
      const double cond = 0.01;

      std::fill_n(bg_conductivities.begin(),nbglayers,cond);
      std::fill_n(bg_thicknesses.begin(),nbglayers,500.0);
      bg_conductivities.back() = 0.1;

      Model.SetBackgroundConductivities(bg_conductivities);
      Model.SetBackgroundThicknesses(bg_thicknesses);
      Model.SetFrequencies().push_back(freqhigh);
      Model.SetFrequencies().push_back(freqlow);
      const size_t nfreq = Model.GetFrequencies().size();

      jiba::OneDMTCalculator Calculator;

      jiba::rvec Impedance = Calculator.Calculate(Model);

      const double prec = 1.0;
      std::complex<double> HighZ(Impedance(0),Impedance(1));
      std::complex<double> LowZ(Impedance(2),Impedance(3));
      BOOST_CHECK_CLOSE(1.0/bg_conductivities.front(),jiba::AppRes(HighZ,freqhigh),prec);
      BOOST_CHECK_CLOSE(1.0/bg_conductivities.back(),jiba::AppRes(LowZ,freqlow),prec);
      BOOST_CHECK_CLOSE(45.0,jiba::ImpedancePhase(HighZ),prec);
      BOOST_CHECK_CLOSE(45.0,jiba::ImpedancePhase(LowZ),prec);
    }

  BOOST_AUTO_TEST_CASE (OneDMT_basic_deriv_test)
    {
      const size_t nbglayers = 7;
      const size_t nfreq = 11;
      jiba::X3DModel Model;
      std::vector<double> bg_thicknesses(nbglayers),bg_conductivities(nbglayers);

      const double freqhigh = 100.0;
      const double freqlow = 0.00001;
      const double cond = 0.01;

      std::fill_n(bg_conductivities.begin(),nbglayers,cond);
      std::fill_n(bg_thicknesses.begin(),nbglayers,500.0);

      Model.SetBackgroundConductivities(bg_conductivities);
      Model.SetBackgroundThicknesses(bg_thicknesses);
      for (size_t i = 0; i< nfreq; ++i)
        {
          double logfreq = std::log10(freqlow)+(std::log10(freqhigh) - std::log10(freqlow))/nfreq * i;
          Model.SetFrequencies().push_back(std::pow(10.0,logfreq));
        }



      jiba::X3DModel TrueModel(Model);

      for (size_t i = 0; i < nbglayers; ++i)
        {
          bg_conductivities.at(i) = std::exp(double(i)-double(nbglayers));
        }

      TrueModel.SetBackgroundConductivities(bg_conductivities);

      jiba::OneDMTCalculator Calculator;
      jiba::rvec Observed = Calculator.Calculate(TrueModel);

      jiba::OneDMTObjective Objective;
      Objective.SetObservedData(Observed);
      Objective.SetModelGeometry(Model);
      Objective.SetDataError(jiba::ConstructError(Observed, 0.02, 0.0));
      jiba::rvec ModelVec(nbglayers);
      std::copy(Model.GetBackgroundConductivities().begin(),Model.GetBackgroundConductivities().end(),ModelVec.begin());
      double misfit = Objective.CalcMisfit(ModelVec);
      BOOST_CHECK(misfit > 0.0);
      jiba::rvec Gradient = Objective.CalcGradient(ModelVec);

      std::ofstream outfile("grad.comp");

      for (size_t index = 0; index < nbglayers; ++index)
        {
          double delta = ModelVec(index) * 0.001;
          jiba::rvec Forward(ModelVec);
          jiba::rvec Backward(ModelVec);
          Forward(index) += delta;
          Backward(index) -= delta;
          double ForFDGrad = (Objective.CalcMisfit(Forward) - misfit)/(delta);
          double BackFDGrad = (misfit - Objective.CalcMisfit(Backward))/delta;
          double CentFDGrad = (ForFDGrad + BackFDGrad)/2.0;
          BOOST_CHECK(Between(ForFDGrad,BackFDGrad,Gradient(index)) || (CentFDGrad - Gradient(index))/CentFDGrad < 0.01);
          outfile << index << " " << bg_conductivities.at(index) << " " << ModelVec(index) << " " << ForFDGrad << " "<< BackFDGrad << " " << (ForFDGrad + BackFDGrad)/2.0 << " " << Gradient(index) << std::endl;
        }
    }
  BOOST_AUTO_TEST_SUITE_END()
