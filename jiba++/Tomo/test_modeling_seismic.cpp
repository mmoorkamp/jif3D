#define BOOST_TEST_MODULE SeismicModel test
#define BOOST_TEST_MAIN ...
#include <boost/test/included/unit_test.hpp>
#include <stdlib.h>
#include <time.h>
#include <vector>
#include "Podvin.h"
#include "modeling_seismic.h"
#include <boost/test/floating_point_comparison.hpp>

BOOST_AUTO_TEST_SUITE( Seismic_Test_Suite )

BOOST_AUTO_TEST_CASE (memory_test)
    {
      // check 0 length allocation
      BOOST_CHECK(jiba::memory(NULL,0,1,"test") == NULL);
      // check simple allocation
      size_t length = 100;
      char * mem_test = jiba::memory(NULL,length, 1, "test");
      BOOST_CHECK( mem_test != NULL);
      for (size_t i = 0; i < length; ++i)
        mem_test[i] = i;
      // check reallocation
      mem_test = jiba::memory(mem_test, 2 * length, 1, "test");
      for (size_t i = 0; i < length; ++i)
        {
          BOOST_CHECK( mem_test[i] == i);
        }
      std::vector<char> dummy;
      size_t failsize = dummy.max_size() * 2; // this should be enough to fail
      BOOST_CHECK_THROW(jiba::memory(NULL,failsize,1,"throw"),std::runtime_error);
      BOOST_CHECK_THROW(jiba::memory(mem_test,failsize,1,"throw"),std::runtime_error);
    }

  BOOST_AUTO_TEST_CASE(podvin3d_test)
    {
      const int nx = 2;
      const int ny = 2;
      const int nz = 2;
      const int nparam = nx * ny * nz;
      float HS[nparam];
      float T[nparam];
      float XS = 0.5;
      float YS = 0.5;
      float ZS = 0.5;
      float HS_EPS_INIT = 0.001;
      int MSG = 2;
      std::fill_n(HS, nparam, 1.0);
      std::fill_n(T, nparam, 0.0);
      int status = time_3d(HS, T, nx, ny, nz, XS, YS, ZS, HS_EPS_INIT, MSG);

      std::cout << T[0] << " " << sqrt(XS * XS + YS * YS + ZS * ZS)
          << std::endl;
    }

  BOOST_AUTO_TEST_CASE(basic_forward_ray_test)
    {

      jiba::GEOMETRY geo;
      jiba::GRID_STRUCT grid;
      jiba::DATA_STRUCT data;
      jiba::RP_STRUCT raypath;

      const size_t ncells = 3;
      grid.nx = ncells;
      grid.ny = ncells;
      grid.nz = ncells;
      grid.h = 70;
      grid.nborder = 0;
      const size_t ngrid = pow(ncells, 3);
      std::fill_n(grid.org, 3, 0.0);
      grid.slow = new double[ngrid];
      std::fill_n(grid.slow, ngrid, 1.0);
      grid.border_index = new int[ngrid];
      std::fill_n(grid.border_index, ngrid, 0);

      data.ndata_seis = 1;
      data.ndata_seis_act = 1;
      data.sno = new int[1];
      data.rno = new int[1];
      data.sno[0] = 0;
      data.rno[0] = 1;
      data.xdist = new float[1];
      data.ydist = new float[1];
      data.zdist = new float[1];
      data.xdist[0] = 50;
      data.ydist[0] = 50;
      data.zdist[0] = 50;
      data.shots = new int[1];
      data.recs = new int[1];
      data.lshots = new int[1];
      data.lrecs = new int[1];

      data.shots[0] = 0;
      data.recs[0] = 1;
      data.lshots[0] = 1;
      data.lrecs[0] = 1;
      data.tcalc = new double[1];

      geo.nrec = 1;
      geo.nshot = 1;
      geo.x = new float[geo.nrec + geo.nshot];
      geo.y = new float[geo.nrec + geo.nshot];
      geo.z = new float[geo.nrec + geo.nshot];
      geo.x[0] = 10.0;
      geo.y[0] = 10.0;
      geo.z[0] = 0.0;
      geo.x[1] = 60.0;
      geo.y[1] = 60.0;
      geo.z[1] = 50.0;
      double dist = sqrt(pow(geo.x[0] - geo.x[1],2) + pow(geo.y[0] - geo.y[1],2) + pow(geo.z[0] - geo.z[1],2));
      ForwardModRay(geo, grid, &data, &raypath, 0);
      std::cout << data.tcalc[0] << " " << dist << std::endl;
      double dummy = 0.0;
    }
BOOST_AUTO_TEST_SUITE_END()
