#define BOOST_TEST_MODULE SeismicModel test
#define BOOST_TEST_MAIN ...
#include <boost/test/included/unit_test.hpp>
#include <stdlib.h>
#include <time.h>
#include <vector>
#include "Podvin.h"
#include "modeling_seismic.h"
#include "TomographyCalculator.h"
#include <boost/test/floating_point_comparison.hpp>
#include <boost/random/linear_congruential.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>

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

  BOOST_AUTO_TEST_CASE(interpolate_test)
    {

      const size_t ncells = 5;
      const size_t npos = 10;
      boost::minstd_rand generator(42u);
      boost::uniform_real<> uni_dist(0, ncells);
      boost::variate_generator<boost::minstd_rand&, boost::uniform_real<> >
          Pos(generator, uni_dist);

      jiba::GRID_STRUCT grid;
      grid.nx = ncells;
      grid.ny = ncells;
      grid.nz = ncells;
      grid.h = 5;
      grid.nborder = 0;
      const size_t ngrid = pow(ncells + 1, 3);
      std::fill_n(grid.org, 3, grid.h / 2.0);
      grid.slow = new double[ngrid];
      std::fill_n(grid.slow, ngrid, 1.0 * grid.h);
      float *data = new float[ngrid];
      float xcomp = Pos();
      float ycomp = Pos();
      float zcomp = Pos();
      for (size_t i = 0; i < ncells + 1; ++i)
        for (size_t j = 0; j < ncells + 1; ++j)
          for (size_t k = 0; k < ncells + 1; ++k)
            {
              data[i * (ncells + 1) * (ncells + 1) + j * (ncells + 1) + k]
                  = xcomp * i + ycomp * j + zcomp * k;
            }
      for (size_t i = 0; i < npos; ++i)
        {
          float xpos = Pos();
          float ypos = Pos();
          float zpos = Pos();
          float inter = jiba::interpolate(xpos, ypos, zpos, grid, data);
          float exact = xcomp * xpos + ycomp * ypos + zcomp * zpos;
          BOOST_CHECK_CLOSE(inter,exact,0.001);
        }
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
      jiba::RP_STRUCT *raypath;

      const size_t ncells = 5;
      grid.nx = ncells;
      grid.ny = ncells;
      grid.nz = ncells;
      grid.h = 35;
      grid.nborder = 0;
      const size_t ngrid = pow(ncells + 1, 3);
      std::fill_n(grid.org, 3, grid.h / 2.0);
      grid.slow = new double[ngrid];
      std::fill_n(grid.slow, ngrid, 1.0 * grid.h);

      data.ndata_seis = 1;
      data.ndata_seis_act = 1;
      data.sno = new int[1];
      data.rno = new int[1];
      data.sno[0] = 1;
      data.rno[0] = 2;
      data.tcalc = new double[1];

      geo.nrec = 1;
      geo.nshot = 1;
      geo.x = new float[geo.nrec + geo.nshot];
      geo.y = new float[geo.nrec + geo.nshot];
      geo.z = new float[geo.nrec + geo.nshot];
      geo.x[0] = 75.0;
      geo.y[0] = 100.0;
      geo.z[0] = 100.0;
      geo.x[1] = 157.0;
      geo.y[1] = 157.0;
      geo.z[1] = 157.0;
      double dist = sqrt(pow(geo.x[0] - geo.x[1], 2) + pow(geo.y[0] - geo.y[1],
          2) + pow(geo.z[0] - geo.z[1], 2));

      raypath = new jiba::RP_STRUCT[1];

      ForwardModRay(geo, grid, &data, raypath, 0);
      std::cout << "Time: " << data.tcalc[0] / 1000.0 << " Dist: " << dist
          << " Rel. Error: " << (data.tcalc[0] / 1000.0 - dist) / dist
          << std::endl;
      double totallength = 0.0;
      for (size_t i = 0; i < raypath[0].nray; ++i)
        {
          std::cout << "Ray x: " << raypath[0].x[i] << " Ray y: "
              << raypath[0].y[i] << " Ray z: " << raypath[0].z[i];
          std::cout << " Length: " << raypath[0].len[i] << " Position: "
              << raypath[0].ele[i] << std::endl;
          totallength += raypath[0].len[i];
        }
      std::cout << "Total length: " << totallength << " Raypath time: " << totallength * grid.h << std::endl;

      jiba::TomographyCalculator Calculator;
      jiba::ThreeDSeismicModel Model;
      Model.SetCellSize(grid.h,ncells,ncells,ncells);
      Model.AddMeasurementPoint(157.0,157.0,157.0);
      Model.AddSource(75.0,100.0,100.0);
      std::fill_n(Model.SetSlownesses().origin(),pow(ncells,3),1.0);
      Model.AddMeasurementConfiguration(0,0);
      jiba::rvec time(Calculator.Calculate(Model));
      std::cout << "Bjoern: " << data.tcalc[0] << " C++: " << time(0) << std::endl;
    }
BOOST_AUTO_TEST_SUITE_END()
