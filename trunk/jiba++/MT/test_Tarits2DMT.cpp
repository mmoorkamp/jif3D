#define BOOST_TEST_MODULE Tarits2DModel test
#define BOOST_TEST_MAIN ...
#include <boost/test/included/unit_test.hpp>
#include <stdlib.h>
#include <time.h>
#include <vector>
#include <Tarits2DMT.h>
#include <boost/test/floating_point_comparison.hpp>

BOOST_AUTO_TEST_SUITE( MT2D_Test_Suite )

BOOST_AUTO_TEST_CASE(HPOL_test)
  {
    double T = 1.0;
    int nx = 10;
    int nz = 10;
    int nionos = 1;
    int natmos = 1;
    double XSizes[nx], ZSizes[nz];
    int modelsize = nx*(nz-nionos-natmos); 
    double Rho[modelsize], Hx_real[modelsize], Hx_imag[modelsize], Ey_real[modelsize], Ey_imag[modelsize];
    double rionos = 1.0;
    epol_(&T, &nx, &nz, &nionos, &natmos, XSizes, ZSizes, Rho, &rionos, Hx_real,Hx_imag,Ey_real,Ey_imag); 
  }

BOOST_AUTO_TEST_SUITE_END()
