//============================================================================
// Name        : test_ReadWriteX3D.cpp
// Author      : Feb 17, 2009
// Version     :
// Copyright   : 2009, mmoorkamp
//============================================================================

#define BOOST_TEST_MODULE OneDMTCalculator test
#define BOOST_TEST_MAIN ...
#include "../Global/Jif3DTesting.h"
#include <boost/test/included/unit_test.hpp>
#include <boost/test/tools/floating_point_comparison.hpp>

#include "../Global/Noise.h"
#include "../Global/NumUtil.h"
#include "OneDMTCalculator.h"
#include "OneDMTObjective.h"
#include "MTEquations.h"
#include "X3DModel.h"

BOOST_AUTO_TEST_SUITE (OneDMTCalculator_Suite)

bool Between(const double limit1, const double limit2, const double value)
  {
    const double upper = std::max(limit1, limit2);
    const double lower = std::min(limit1, limit2);
    return (lower <= value) && (upper >= value);
  }

BOOST_AUTO_TEST_CASE (OneD_forward_hs_test)
  {
    const size_t nbglayers = 5;
    const size_t nfreq = 6;
    jif3D::X3DModel Model;
    std::vector<double> bg_thicknesses(nbglayers), bg_conductivities(nbglayers);
    jif3D::MTData Data;
    const double freq = 1.0;
    const double cond = 0.01;

    std::fill_n(bg_conductivities.begin(),nbglayers,cond);
    std::fill_n(bg_thicknesses.begin(),nbglayers,200.0);

    Model.SetBackgroundConductivities(bg_conductivities);
    Model.SetBackgroundThicknesses(bg_thicknesses);
    std::vector<double> Freqs;
    for (size_t i = 0; i < nfreq; ++i)
      {
        Freqs.push_back(freq*(i+1));
      }
    Data.SetFrequencies(Freqs);

    jif3D::OneDMTCalculator Calculator;

    std::vector<double> Impedance = Calculator.Calculate(Model, Data);

    const double prec = 0.05;
    for (size_t i = 0; i < nfreq; ++i)
      {
        std::complex<double> HsImp = jif3D::ImpedanceHalfspace(freq*(i+1),cond);
        BOOST_CHECK_CLOSE(Impedance.at(i*2),HsImp.real(),prec);
        BOOST_CHECK_CLOSE(Impedance.at(i*2+1),HsImp.imag(),prec);
      }

  }

BOOST_AUTO_TEST_CASE (OneD_twolayer_test)
  {
    const size_t nbglayers = 5;

    jif3D::X3DModel Model;
    jif3D::MTData Data;

    std::vector<double> bg_thicknesses(nbglayers), bg_conductivities(nbglayers);

    const double freqhigh = 100.0;
    const double freqlow = 0.00001;
    const double cond = 0.01;

    std::fill_n(bg_conductivities.begin(),nbglayers,cond);
    std::fill_n(bg_thicknesses.begin(),nbglayers,500.0);
    bg_conductivities.back() = 0.1;

    Model.SetBackgroundConductivities(bg_conductivities);
    Model.SetBackgroundThicknesses(bg_thicknesses);
    Data.SetFrequencies(
          { freqhigh,freqlow});

    jif3D::OneDMTCalculator Calculator;

    std::vector<double> Impedance = Calculator.Calculate(Model, Data);

    const double prec = 1.0;
    std::complex<double> HighZ(Impedance.at(0), Impedance.at(1));
    std::complex<double> LowZ(Impedance.at(2), Impedance.at(3));
    BOOST_CHECK_CLOSE(1.0/bg_conductivities.front(),jif3D::AppRes(HighZ,freqhigh),prec);
    BOOST_CHECK_CLOSE(1.0/bg_conductivities.back(),jif3D::AppRes(LowZ,freqlow),prec);
    BOOST_CHECK_CLOSE(45.0,jif3D::ImpedancePhase(HighZ),prec);
    BOOST_CHECK_CLOSE(45.0,jif3D::ImpedancePhase(LowZ),prec);
  }

BOOST_AUTO_TEST_CASE (OneDMT_basic_deriv_test)
  {
    const size_t nbglayers = 7;
    const size_t nfreq = 11;
    jif3D::X3DModel Model;
    jif3D::MTData Data;

    const double freqhigh = 100.0;
    const double freqlow = 0.00001;
    const double cond = 0.01;

    std::vector<double> bg_thicknesses(nbglayers, 500.0), bg_conductivities(nbglayers,
        cond);
    Model.SetBackgroundConductivities(bg_conductivities);
    Model.SetBackgroundThicknesses(bg_thicknesses);
    std::vector<double> Freqs;
    for (size_t i = 0; i< nfreq; ++i)
      {
        double logfreq = std::log10(freqlow)+(std::log10(freqhigh) - std::log10(freqlow))/nfreq * i;
        Freqs.push_back(std::pow(10.0,logfreq));
      }
    Data.SetFrequencies(Freqs);

    jif3D::X3DModel TrueModel(Model);

    for (size_t i = 0; i < nbglayers; ++i)
      {
        bg_conductivities.at(i) = std::exp(double(i)-double(nbglayers));
      }

    TrueModel.SetBackgroundConductivities(bg_conductivities);

    jif3D::OneDMTCalculator Calculator;
    std::vector<double> Observed = Calculator.Calculate(TrueModel, Data);
    std::vector<double> Error(Observed.size(), 0.0);
    Data.SetDataAndErrors(std::vector<double>(Observed.begin(),Observed.end()),Error);
    jif3D::OneDMTObjective Objective;
    Objective.SetObservedData(Data);
    Objective.SetModelGeometry(Model);
    Objective.SetDataError(jif3D::ConstructError(std::vector<double>(Observed.begin(),Observed.end()), Error, 0.02, 0.0));
    jif3D::rvec ModelVec(nbglayers);
    std::copy(Model.GetBackgroundConductivities().begin(),Model.GetBackgroundConductivities().end(),ModelVec.begin());
    double misfit = Objective.CalcMisfit(ModelVec);
    BOOST_CHECK(misfit > 0.0);
    jif3D::rvec Gradient = Objective.CalcGradient(ModelVec);

    std::ofstream outfile("grad1d.comp");

    for (size_t index = 0; index < nbglayers; ++index)
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
        outfile << index << ForFDGrad << " "<< BackFDGrad << " " << (ForFDGrad + BackFDGrad)/2.0 << " " << Gradient(index) << std::endl;
      }
  }
BOOST_AUTO_TEST_SUITE_END()
