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

#include "../Global/NumUtil.h"
#include "X3DMTCalculator.h"
#include "X3DModel.h"
#include "MTEquations.h"
#include "MT2DForward.h"
#include "ReadWriteImpedances.h"

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
      Model.SetHorizontalCellSize(deltax,deltay,xsize,ysize);

      std::fill_n(Model.SetZCellSizes().origin(),zsize,deltaz);
      std::fill_n(Model.SetConductivities().origin(),xsize*ysize*zsize,0.01);
      std::fill_n(bg_conductivities.begin(),nbglayers,0.01);
      std::fill_n(bg_thicknesses.begin(),nbglayers,100.0);

      Model.SetBackgroundConductivities(bg_conductivities);
      Model.SetBackgroundThicknesses(bg_thicknesses);
      Model.SetFrequencies().push_back(10.0);

      jiba::X3DMTCalculator Calculator;
      jiba::rvec Impedance = Calculator.Calculate(Model);

      std::complex<double> HsImp = jiba::ImpedanceHalfspace(10.0,0.01);
      for (size_t i = 0; i < xsize * ysize; ++i)
        {
          BOOST_CHECK_CLOSE(Impedance(i*8+2),HsImp.real(),0.01);
          BOOST_CHECK_CLOSE(Impedance(i*8+3),HsImp.imag(),0.01);
          BOOST_CHECK_CLOSE(Impedance(i*8+4),-HsImp.real(),0.01);
          BOOST_CHECK_CLOSE(Impedance(i*8+5),-HsImp.imag(),0.01);
        }
    }
  BOOST_AUTO_TEST_CASE (X3D_forward_2D_test)
    {
      const size_t xsize = 100;
      const size_t ysize = 100;
      const size_t zsize = 15;
      const size_t nbglayers = 15;
      jiba::X3DModel Model;

      Model.SetConductivities().resize(boost::extents[xsize][ysize][zsize]);
      std::vector<double> bg_thicknesses(nbglayers),bg_conductivities(nbglayers);

      const double deltax = 100.0;
      const double deltay = 100.0;
      const double deltaz = 100.0;
      const double Period = 1.0;
      const double bg_cond = 0.1;
      const double anom_cond = 1.0;

      Model.SetHorizontalCellSize(deltax,deltay,xsize,ysize);
      Model.SetZCellSizes().resize(boost::extents[zsize]);

      double currsize = 200.0;
      for (int i = 0; i < zsize; ++i)
        {
          bg_thicknesses[i] = floor(currsize);
          Model.SetZCellSizes()[i] = floor(currsize);
          currsize *= 1.1;
        }
      //std::fill_n(bg_thicknesses.begin(),nbglayers,deltaz);
      //std::fill_n(Model.SetZCellSizes().origin(),zsize,deltaz);

      std::fill_n(Model.SetConductivities().origin(),xsize*ysize*zsize,bg_cond);
      std::fill_n(bg_conductivities.begin(),nbglayers,bg_cond * 1.001);

      for (size_t i = 40; i < 60; ++i )
      for (size_t j = 0; j < ysize; ++j)
      for (size_t k = 0; k < zsize; ++k)
        {
          Model.SetConductivities()[i][j][k] = anom_cond;
        }

      Model.SetBackgroundConductivities(bg_conductivities);
      Model.SetBackgroundThicknesses(bg_thicknesses);
      Model.SetFrequencies().push_back(1.0/Period);

      jiba::X3DMTCalculator Calculator;
      jiba::rvec Impedance3D = Calculator.Calculate(Model);
      jiba::rvec Freq(Model.GetFrequencies().size());
      std::copy(Model.GetFrequencies().begin(),Model.GetFrequencies().end(),Freq.begin());
      std::vector<double> XCoord(xsize*ysize),YCoord(xsize*ysize),ZCoord(xsize*ysize);
      std::fill_n(ZCoord.begin(),xsize*ysize,0.0);
      for (size_t i = 0; i < Model.GetXCoordinates().size(); ++i)
        {
          for (size_t j = 0; j < Model.GetYCoordinates().size(); ++j)
            {
              XCoord[i] = Model.GetXCoordinates()[i];
              YCoord[i] = Model.GetYCoordinates()[i];

            }
        }
      jiba::WriteImpedancesToNetCDF("imp2D.nc",Freq,XCoord,YCoord,ZCoord,Impedance3D);
      jiba::rvec Imp3DProfile(xsize*8);
      for (size_t i = 0; i < xsize; ++i)
        {
          std::copy(Impedance3D.begin()+i*ysize*8+400,Impedance3D.begin()+i*ysize*8+8+400,Imp3DProfile.begin()+i*8);
        }
      std::vector<double> XCoordProf(xsize),YCoordProf(xsize),ZCoordProf(xsize);
      std::fill_n(ZCoordProf.begin(),xsize,0.0);
      std::fill_n(YCoordProf.begin(),xsize,0.0);
      std::generate_n(XCoordProf.begin(),xsize,jiba::IntSequence(0));
      jiba::WriteImpedancesToNetCDF("imp3Dprof.nc",Freq,XCoordProf,YCoordProf,ZCoordProf,Imp3DProfile);
      Model.WriteVTK("mod2D.vtk");

      jiba::MT2DForward Forward2D;
      const size_t nx2D = 100;
      const size_t nz2D = 50;
      Forward2D.SetXSizes().resize(boost::extents[nx2D]);
      Forward2D.SetZSizes().resize(boost::extents[nz2D]);
      currsize = 10.0;
      std::fill_n(Forward2D.SetXSizes().origin(),nx2D,100.0);
      for (int i = 0; i < nz2D; ++i)
        {
          Forward2D.SetZSizes()[i] = currsize;
          currsize *= 1.1;
        }
      Forward2D.SetResistivities().resize(boost::extents[nx2D][nz2D]);
      std::fill_n(Forward2D.SetResistivities().origin(),nx2D*nz2D,1./bg_cond);
      for (size_t i = 40; i < 60; ++i )
      for (size_t k = 0; k < nz2D; ++k)
      Forward2D.SetResistivities()[i][k] = 1./anom_cond;

      std::vector<double> Periods(1,Period);
      Forward2D.CalcEpol(Periods);
      Forward2D.CalcBpol(Periods);

      for (size_t i = 0; i < nx2D; ++i)
        {
          std::complex<double> Hx(Forward2D.GetHx_real()[0][i][0],Forward2D.GetHx_imag()[0][i][0]);
          std::complex<double> Ey(Forward2D.GetEy_real()[0][i][0],Forward2D.GetEy_imag()[0][i][0]);
          std::complex<double> Zyx(Ey/Hx);

          std::complex<double> Hy(Forward2D.GetHy_real()[0][i][0],Forward2D.GetHy_imag()[0][i][0]);
          std::complex<double> Ex(Forward2D.GetEx_real()[0][i][0],Forward2D.GetEx_imag()[0][i][0]);
          std::complex<double> Zxy(Ex/Hy);
          BOOST_CHECK_CLOSE(Zxy.real(),Imp3DProfile(i*8+2),1.0);
          BOOST_CHECK_CLOSE(Zxy.imag(),Imp3DProfile(i*8+3),1.0);
          BOOST_CHECK_CLOSE(Zyx.real(),Imp3DProfile(i*8+4),1.0);
          BOOST_CHECK_CLOSE(Zyx.imag(),Imp3DProfile(i*8+5),1.0);
          Imp3DProfile(i*8+2) = Zxy.real();
          Imp3DProfile(i*8+3) = Zxy.imag();
          Imp3DProfile(i*8+4) = Zyx.real();
          Imp3DProfile(i*8+5) = Zyx.imag();
          BOOST_CHECK(fabs(Imp3DProfile(i*8+0)/Imp3DProfile(i*8+2)) < 1e-3);
          BOOST_CHECK(fabs(Imp3DProfile(i*8+1)/Imp3DProfile(i*8+3)) < 1e-3);
          BOOST_CHECK(fabs(Imp3DProfile(i*8+6)/Imp3DProfile(i*8+4)) < 1e-3);
          BOOST_CHECK(fabs(Imp3DProfile(i*8+7)/Imp3DProfile(i*8+5)) < 1e-3);
        }
      jiba::WriteImpedancesToNetCDF("imp2Dprof.nc",Freq,XCoordProf,YCoordProf,ZCoordProf,Imp3DProfile);
    }

  BOOST_AUTO_TEST_SUITE_END()
