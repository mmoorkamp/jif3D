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
#include "X3DMTCalculator.h"
#include "X3DModel.h"
#include "OneDMTCalculator.h"
#include "MTEquations.h"
#include "MT2DForward.h"
#include "ReadWriteImpedances.h"
#include "MTData.h"

BOOST_AUTO_TEST_SUITE (X3DCalculator_Suite)

BOOST_AUTO_TEST_CASE (X3D_forward_hs_test)
  {
    //create a 3D version of a 1D model
    const size_t xsize = 20;
    const size_t ysize = 20;
    const size_t zsize = 10;
    const size_t nbglayers = 5;
    jif3D::X3DModel Model;
    jif3D::MTData Data;

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

    jif3D::X3DMTCalculator Calculator;
    for (size_t i = 1; i < xsize / 2; ++i)
    for (size_t j = 1; j < ysize / 2; ++j)
      {
        Data.AddMeasurementPoint(Model.GetXCoordinates()[i] + deltax / 2.0,
            Model.GetYCoordinates()[j] + deltay / 2.0, 0.0);
      }
    Data.CompleteObject();
    jif3D::rvec Impedance = Calculator.Calculate(Model, Data);
    std::complex<double> HsImp = jif3D::ImpedanceHalfspace(freq, cond);
    const double prec = 0.05;
    for (size_t i = 0; i < Impedance.size() / 8; ++i)
      {
        BOOST_CHECK_CLOSE(Impedance(i * 8 + 2), HsImp.real(), prec);
        BOOST_CHECK_CLOSE(Impedance(i * 8 + 3), HsImp.imag(), prec);
        BOOST_CHECK_CLOSE(Impedance(i * 8 + 4), -HsImp.real(), prec);
        BOOST_CHECK_CLOSE(Impedance(i * 8 + 5), -HsImp.imag(), prec);
      }

  }

BOOST_AUTO_TEST_CASE (X3D_layered_hs_test)
  {
    //create a 3D version of a 1D model
    const size_t xsize = 20;
    const size_t ysize = 20;
    const size_t zsize = 10;
    const size_t nbglayers = 5;
    jif3D::X3DModel Model;
    jif3D::MTData Data;

    Model.SetConductivities().resize(boost::extents[xsize][ysize][zsize]);
    std::vector<double> bg_thicknesses(nbglayers), bg_conductivities(nbglayers);

    const double deltax = 100.0;
    const double deltay = 100.0;
    const double deltaz = 100.0;

    const double cond = 0.01;
    Model.SetHorizontalCellSize(deltax, deltay, xsize, ysize);
    jif3D::ThreeDModelBase::t3DModelDim ZCS(zsize, deltaz);
    Model.SetZCellSizes(ZCS);
    std::fill_n(Model.SetConductivities().origin(), xsize * ysize * zsize, cond);
    std::fill_n(bg_conductivities.begin(), nbglayers, cond);
    std::fill_n(bg_thicknesses.begin(), nbglayers, 200.0);

    bg_conductivities.at(2) = 1.0;
    for (size_t i = 0; i < xsize; ++i)
    for (size_t j = 0; j < ysize; ++j)
      {
        Model.SetConductivities()[i][j][4] = 1.0;
        Model.SetConductivities()[i][j][5] = 1.0;
      }
    Model.SetConductivities()[0][0][9] = 1.01;
    Model.SetBackgroundConductivities(bg_conductivities);
    Model.SetBackgroundThicknesses(bg_thicknesses);
    Data.SetFrequencies(
          { 1000.0, 100.0, 10.0, 1.0, 0.1, 0.01});

    jif3D::X3DMTCalculator Calculator;
    for (size_t i = 1; i < xsize / 2; ++i)
    for (size_t j = 1; j < ysize / 2; ++j)
      {
        Data.AddMeasurementPoint(Model.GetXCoordinates()[i] + deltax / 2.0,
            Model.GetYCoordinates()[j] + deltay / 2.0, 0.0);
      }
    Data.CompleteObject();
    jif3D::rvec Impedance3D = Calculator.Calculate(Model, Data);
    jif3D::rvec Impedance1D = jif3D::OneDMTCalculator().Calculate(Model, Data);
    Data.SetDataAndErrors(std::vector<double>(Impedance3D.begin(),Impedance3D.end()),std::vector<double>(Impedance3D.size(),0.0));
    Data.WriteNetCDF("imp3D.nc");

    const size_t nsites = Data.GetMeasPosX().size();
    const size_t nfreq = Data.GetFrequencies().size();
    const double prec = 0.05;
    for (size_t i = 0; i < nsites; ++i)
      {
        for (size_t j = 0; j < nfreq; ++j)
          {
            BOOST_CHECK_CLOSE(Impedance3D((j * nsites + i) * 8 + 2),
                Impedance1D(2 * j), prec);
            BOOST_CHECK_CLOSE(Impedance3D((j * nsites + i) * 8 + 3),
                Impedance1D(2 * j + 1), prec);
            BOOST_CHECK_CLOSE(Impedance3D((j * nsites + i) * 8 + 4),
                -Impedance1D(2 * j), prec);
            BOOST_CHECK_CLOSE(Impedance3D((j * nsites + i) * 8 + 5),
                -Impedance1D(2 * j + 1), prec);
          }
      }

  }

BOOST_AUTO_TEST_CASE (X3D_forward_2D_test)
  {
    const size_t xsize = 100;
    const size_t ysize = 50;
    const size_t zsize = 20;
    const size_t nbglayers = zsize;
    jif3D::X3DModel Model;
    jif3D::MTData Data;

    Model.SetConductivities().resize(boost::extents[xsize][ysize][zsize]);
    std::vector<double> bg_thicknesses(nbglayers), bg_conductivities(nbglayers);

    const double deltax = 100.0;
    const double deltay = 500.0;
    const double deltaz = 100.0;
    const double Period = 2.0;
    const double bg_cond = 0.1;
    const double anom_cond = 0.5;

    Model.SetHorizontalCellSize(deltax, deltay, xsize, ysize);
    jif3D::ThreeDModelBase::t3DModelDim ZCS(zsize, deltaz);

    double currsize = deltaz;
    for (size_t i = 0; i < zsize; ++i)
      {
        bg_thicknesses[i] = floor(currsize);
        ZCS[i] = floor(currsize);
        currsize *= 1.1;
      }
    Model.SetZCellSizes(ZCS);
    //std::fill_n(bg_thicknesses.begin(), nbglayers, bg_thick);
    //std::fill_n(Model.SetZCellSizes().origin(),zsize,deltaz);

    std::fill_n(Model.SetConductivities().origin(), xsize * ysize * zsize, bg_cond);
    std::fill_n(bg_conductivities.begin(), nbglayers, bg_cond);

    for (size_t i = 40; i < 60; ++i)
    for (size_t j = 0; j < ysize; ++j)
    for (size_t k = 0; k < zsize; ++k)
      {
        Model.SetConductivities()[i][j][k] = anom_cond;
      }

    Model.SetBackgroundConductivities(bg_conductivities);
    Model.SetBackgroundThicknesses(bg_thicknesses);
    Data.SetFrequencies(
          { 1.0 / Period});
    for (size_t i = 1; i < xsize - 1; ++i)
      {
        Data.AddMeasurementPoint(Model.GetXCoordinates()[i] + deltax / 2.0,
            (ysize * deltay) / 2.0 + deltay / 2.0, 0.0);
      }
    Data.CompleteObject();

    jif3D::X3DMTCalculator Calculator;
    std::cout << "Calculating 3D forward " << std::endl;
    jif3D::rvec Impedance3D = Calculator.Calculate(Model, Data);

    Data.SetDataAndErrors(std::vector<double>(Impedance3D.begin(),Impedance3D.end()),std::vector<double>(Impedance3D.size(),0.0));
    Data.WriteNetCDF("imp3Dprof.nc");

    jif3D::rvec Imp3DProfile(Impedance3D);

    Model.WriteVTK("mod2D.vtk");

    jif3D::MT2DForward Forward2D;
    const size_t nx2D = 100;
    const size_t nz2D = 50;
    Forward2D.SetXSizes().resize(boost::extents[nx2D]);
    Forward2D.SetZSizes().resize(boost::extents[nz2D]);
    double currzsize = 10.0;
    double currxsize = 100.0;
    std::fill_n(Forward2D.SetXSizes().origin(), nx2D, 100.0);
    for (size_t i = 10; i > 0; --i)
      {
        Forward2D.SetXSizes()[i] = currxsize;
        currxsize *= 1.1;
      }
    for (size_t i = 90; i < nx2D; ++i)
      {
        Forward2D.SetXSizes()[i] = currxsize;
        currxsize *= 1.1;
      }
    for (size_t i = 0; i < nz2D; ++i)
      {
        Forward2D.SetZSizes()[i] = currzsize;
        currzsize *= 1.1;
      }
    Forward2D.SetResistivities().resize(boost::extents[nx2D][nz2D]);
    std::fill_n(Forward2D.SetResistivities().origin(), nx2D * nz2D, 1. / bg_cond);
    for (size_t i = 40; i < 60; ++i)
    for (size_t k = 0; k < nz2D; ++k)
    Forward2D.SetResistivities()[i][k] = 1. / anom_cond;

    std::vector<double> Periods(1, Period);
    std::cout << "Calculating 2D forward " << std::endl;
    Forward2D.CalcEpol(Periods);
    Forward2D.CalcBpol(Periods);
    std::ofstream diffile("2dmtdiff.out");
    double zxyrchi = 0, zxyichi = 0, zyxrchi = 0, zyxichi = 0;
    for (size_t i = 1; i < nx2D - 1; ++i)
      {
        std::complex<double> Hx(Forward2D.GetHx_real()[0][i][0],
            Forward2D.GetHx_imag()[0][i][0]);
        std::complex<double> Ey(Forward2D.GetEy_real()[0][i][0],
            Forward2D.GetEy_imag()[0][i][0]);
        std::complex<double> Zyx(Ey / Hx);

        std::complex<double> Hy(Forward2D.GetHy_real()[0][i][0],
            Forward2D.GetHy_imag()[0][i][0]);
        std::complex<double> Ex(Forward2D.GetEx_real()[0][i][0],
            Forward2D.GetEx_imag()[0][i][0]);
        std::complex<double> Zxy(Ex / Hy);

        zxyrchi += pow((Zxy.real() - Imp3DProfile((i - 1) * 8 + 2)) / Zxy.real(), 2);
        zxyichi += pow((Zxy.imag() - Imp3DProfile((i - 1) * 8 + 3)) / Zxy.imag(), 2);
        zyxrchi += pow((Zyx.real() - Imp3DProfile((i - 1) * 8 + 4)) / Zyx.real(), 2);
        zyxichi += pow((Zyx.imag() - Imp3DProfile((i - 1) * 8 + 5)) / Zyx.imag(), 2);
        diffile << i << " "
        << (Zxy.real() - Imp3DProfile((i - 1) * 8 + 2)) / Zxy.real() << " "
        << (Zxy.imag() - Imp3DProfile((i - 1) * 8 + 3)) / Zxy.imag() << " "
        << (Zyx.real() - Imp3DProfile((i - 1) * 8 + 4)) / Zyx.real() << " "
        << (Zyx.imag() - Imp3DProfile((i - 1) * 8 + 5)) / Zyx.imag() << "\n";
        //6 percent is not that close, but we have some differences near the boundary
        //of the anomaly
        BOOST_CHECK_CLOSE(Zxy.real(), Imp3DProfile((i - 1) * 8 + 2), 6.0);
        BOOST_CHECK_CLOSE(Zxy.imag(), Imp3DProfile((i - 1) * 8 + 3), 6.0);
        BOOST_CHECK_CLOSE(Zyx.real(), Imp3DProfile((i - 1) * 8 + 4), 6.0);
        BOOST_CHECK_CLOSE(Zyx.imag(), Imp3DProfile((i - 1) * 8 + 5), 6.0);

        BOOST_CHECK(fabs(Imp3DProfile((i - 1) * 8 + 0) / Zxy.real()) < 1e-3);
        BOOST_CHECK(fabs(Imp3DProfile((i - 1) * 8 + 1) / Zxy.imag()) < 1e-3);
        BOOST_CHECK(fabs(Imp3DProfile((i - 1) * 8 + 6) / Zyx.real()) < 1e-3);
        BOOST_CHECK(fabs(Imp3DProfile((i - 1) * 8 + 7) / Zyx.imag()) < 1e-3);
      }
    BOOST_CHECK(sqrt(zxyrchi / nx2D) < 0.02);
    BOOST_CHECK(sqrt(zxyichi / nx2D) < 0.02);
    BOOST_CHECK(sqrt(zyxrchi / nx2D) < 0.02);
    BOOST_CHECK(sqrt(zyxichi / nx2D) < 0.02);

    std::cout << sqrt(zxyrchi / nx2D) << " " << sqrt(zxyichi / nx2D) << " "
    << sqrt(zyxrchi / nx2D) << " " << sqrt(zyxichi / nx2D) << "\n";
    Data.SetDataAndErrors(std::vector<double>(Imp3DProfile.begin(),Imp3DProfile.end()),std::vector<double>(Impedance3D.size(),0.0));
    Data.WriteNetCDF("imp2Dprof.nc");

  }
BOOST_AUTO_TEST_SUITE_END()
