#define BOOST_TEST_MODULE Tarits2DModel test
#define BOOST_TEST_MAIN ...
#include <boost/test/included/unit_test.hpp>
#include <vector>
#include <complex>
#include "Tarits2DMT.h"
#include "MT2DForward.h"
#include <boost/test/floating_point_comparison.hpp>
#include <boost/lambda/lambda.hpp>

BOOST_AUTO_TEST_SUITE( MT2D_Test_Suite )

BOOST_AUTO_TEST_CASE  (EPOL_test)
    {
      double T = 10.0;
      long nx = 100;
      long nz = 150;
      long nionos = 20;
      long natmos = 20;
      std::vector<double> XSizes(nx,0.0), ZSizes(nz,0.0);
      int modelsize = nx*(nz-nionos-natmos);
      std::vector<double> Rho, Hx_real, Hx_imag, Ey_real, Ey_imag;
      double rionos = 1.0;
      std::fill_n(XSizes.begin(),nx,100.0);
      double currsize = 10.0;
      for (int i = 0; i < nionos+natmos; ++i)
        {
          ZSizes[i] = 100.0;
        }
      for (int i = nionos+natmos; i < nz; ++i)
        {
          ZSizes[i] = currsize;
          currsize *= 1.1;
        }

      std::fill_n(back_inserter(Rho),modelsize,10.0);
      std::fill_n(back_inserter(Hx_real),modelsize,0.0);
      std::fill_n(back_inserter(Hx_imag),modelsize,0.0);
      std::fill_n(back_inserter(Ey_real),modelsize,0.0);
      std::fill_n(back_inserter(Ey_imag),modelsize,0.0);
      epol_(&T, &nx, &nz, &nionos, &natmos, &XSizes[0], &ZSizes[0], &Rho[0], &rionos, &Hx_real[0],&Hx_imag[0],&Ey_real[0],&Ey_imag[0]);

      const size_t index = nx/2 *nz;
      std::complex<double> Hx(Hx_real[index],Hx_imag[index]);
      std::complex<double> Ey(Ey_real[index],Ey_imag[index]);
      std::complex<double> Z(Ey/Hx);
      std::cout << "Z: " << Z << " arg: " << std::arg(Z) * 180.0/M_PI << " rho_a: " << T*abs(Z) * abs(Z)/(8.0e-7*M_PI*M_PI) << std::endl;

    }

  BOOST_AUTO_TEST_CASE(HPOL_test)
    {
      double T = 10.0;
      long nx = 100;
      long nz = 150;
      std::vector<double> XSizes(nx,0.0), ZSizes(nz,0.0);
      int modelsize = nx*nz;
      std::vector<double> Rho, Hy_real, Hy_imag, Ex_real, Ex_imag,Ez_real, Ez_imag;
      std::fill_n(XSizes.begin(),nx,100.0);
      double currsize = 10.0;
      for (int i = 0; i < nz; ++i)
        {
          ZSizes[i] = currsize;
          currsize *= 1.15;
        }

      std::fill_n(back_inserter(Rho),modelsize,10.0);
      std::fill_n(back_inserter(Hy_real),modelsize,0.0);
      std::fill_n(back_inserter(Hy_imag),modelsize,0.0);
      std::fill_n(back_inserter(Ex_real),modelsize,0.0);
      std::fill_n(back_inserter(Ex_imag),modelsize,0.0);
      std::fill_n(back_inserter(Ez_real),modelsize,0.0);
      std::fill_n(back_inserter(Ez_imag),modelsize,0.0);
      hpol_(&T, &nx, &nz, &XSizes[0], &ZSizes[0], &Rho[0], &Hy_real[0],&Hy_imag[0],&Ex_real[0],&Ex_imag[0],&Ez_real[0],&Ez_imag[0]);

      const size_t index = nx/2 *nz;
      std::complex<double> Hy(Hy_real[index],Hy_imag[index]);
      std::complex<double> Ex(Ex_real[index],Ex_imag[index]);
      std::complex<double> Z(Ex/Hy);
      std::cout << "Z: " << Z << " arg: " << std::arg(Z) * 180.0/M_PI << " rho_a: " << T*abs(Z) * abs(Z)/(8.0e-7*M_PI*M_PI) << std::endl;

    }

  BOOST_AUTO_TEST_CASE(MT2DForward_test)
    {
      jif3D::MT2DForward Forward;
      const size_t nx = 100;
      const size_t nz = 100;
      Forward.SetXSizes().resize(boost::extents[nx]);
      Forward.SetZSizes().resize(boost::extents[nz]);
      double currsize = 10.0;
      std::fill_n(Forward.SetXSizes().origin(),nx,100.0);
      for (int i = 0; i < nz; ++i)
        {
          Forward.SetZSizes()[i] = currsize;
          currsize *= 1.1;
        }
      Forward.SetResistivities().resize(boost::extents[nx][nz]);
      std::fill_n(Forward.SetResistivities().origin(),nx*nz,10.0);
      std::vector<double> Periods(2,10.0);
      Periods.at(1) = 100.0;
      Forward.CalcEpol(Periods);
      Forward.CalcBpol(Periods);
      for (int i = 0; i < 2; ++i)
        {
          const size_t index = 50;
          std::complex<double> Hx(Forward.GetHx_real()[i][index][0],Forward.GetHx_imag()[i][index][0]);
          std::complex<double> Ey(Forward.GetEy_real()[i][index][0],Forward.GetEy_imag()[i][index][0]);
          std::complex<double> Zyx(Ey/Hx);
          std::cout << "Z: " << Zyx << " arg: " << std::arg(Zyx) * 180.0/M_PI << " rho_a: " << Periods.at(i)*abs(Zyx) * abs(Zyx)/(8.0e-7*M_PI*M_PI) << std::endl;
          std::complex<double> Hy(Forward.GetHy_real()[i][index][0],Forward.GetHy_imag()[i][index][0]);
          std::complex<double> Ex(Forward.GetEx_real()[i][index][0],Forward.GetEx_imag()[i][index][0]);
          std::complex<double> Zxy(Ex/Hy);
          std::cout << "Z: " << Zxy << " arg: " << std::arg(Zxy) * 180.0/M_PI << " rho_a: " << Periods.at(i)*abs(Zxy) * abs(Zxy)/(8.0e-7*M_PI*M_PI) << std::endl;
        }
    }
  BOOST_AUTO_TEST_SUITE_END()
