
#define BOOST_TEST_MODULE SeismicModel test
#define BOOST_TEST_MAIN ...
#include "../Global/Jif3DTesting.h"
#include <stdlib.h>
#include <time.h>
#include <vector>
#include "../Tomo/PodvinTime3D.h"
#include "../Tomo/modeling_seismic.h"
#include "../Tomo/TomographyCalculator.h"
#include "../Tomo/TomographyData.h"
#include "../Tomo/ReadWriteTomographyData.h"
#include <boost/test/tools/floating_point_comparison.hpp>
#include <boost/random/linear_congruential.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>

BOOST_AUTO_TEST_SUITE (Seismic_Test_Suite)

//test that the interpolation function works for a strictly quadratic function
BOOST_AUTO_TEST_CASE (interpolate_test)
  {

    const size_t ncells = 5;
    const size_t npos = 10;
    boost::minstd_rand generator(42u);
    boost::uniform_real<> uni_dist(0, ncells);
    boost::variate_generator<boost::minstd_rand&, boost::uniform_real<> > Pos(generator,
        uni_dist);

    jif3D::GRID_STRUCT grid;
    grid.nx = ncells;
    grid.ny = ncells;
    grid.nz = ncells;
    grid.h = 5;
    const size_t ngrid = pow(double(ncells + 1), 3);
    //std::fill_n(grid.org, 3, grid.h / 2.0);
    grid.slow.resize(ngrid);
    std::fill_n(grid.slow.begin(), ngrid, 1.0 * grid.h);
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
        float inter = jif3D::interpolate(xpos, ypos, zpos, grid, data);
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
    int MSG = 0;
    std::fill_n(HS, nparam, 1.0);
    std::fill_n(T, nparam, 0.0);
    jif3D::PodvinTime3D()
    .time_3d(HS, T, nx, ny, nz, XS, YS, ZS,
        HS_EPS_INIT, MSG);
    float dist = sqrt(XS * XS + YS * YS + ZS * ZS);
    BOOST_CHECK_CLOSE(T[0],dist,0.001);
  }

BOOST_AUTO_TEST_CASE(basic_forward_ray_test)
  {

    jif3D::GEOMETRY geo;
    jif3D::GRID_STRUCT grid;
    jif3D::DATA_STRUCT data;
    std::vector<jif3D::RP_STRUCT> raypath;

    const size_t ncells = 5;
    grid.nx = ncells;
    grid.ny = ncells;
    grid.nz = ncells;
    grid.h = 35;
    const size_t ngrid = pow(double(ncells + 1), 3);
    const double slow = 0.1;
    grid.slow.resize(ngrid);
    std::fill_n(grid.slow.begin(), ngrid, slow * grid.h);

    data.ndata_seis = 2;
    data.ndata_seis_act = 2;
    data.sno.resize(2);
    data.rno.resize(2);
    data.sno[0] = 1;
    data.rno[0] = 2;
    data.sno[1] = 2;
    data.rno[1] = 1;
    data.tcalc.resize(2);

    geo.nrec = 2;
    geo.nshot = 2;
    geo.x.resize(geo.nrec + geo.nshot);
    geo.y.resize(geo.nrec + geo.nshot);
    geo.z.resize(geo.nrec + geo.nshot);
    geo.x[0] = 75.0;
    geo.y[0] = 100.0;
    geo.z[0] = 40.0;
    geo.x[1] = 157.0;
    geo.y[1] = 157.0;
    geo.z[1] = 157.0;
    double dist = sqrt(
        pow(geo.x[0] - geo.x[1], 2) + pow(geo.y[0] - geo.y[1], 2)
            + pow(geo.z[0] - geo.z[1], 2));

    raypath.resize(2);

    ForwardModRay(geo, grid, data, raypath);
    double relerror = (data.tcalc[0] - dist * slow) / (dist * slow);
    BOOST_CHECK(std::abs(relerror) < 0.05 );

    double totallength = 0.0;
    for (size_t i = 0; i < raypath[0].nray + 1; ++i)
      {
        totallength += raypath[0].len[i];
      }
    BOOST_CHECK(std::abs(dist-totallength * grid.h) < slow * grid.h);
    jif3D::TomographyCalculator Calculator;
    jif3D::ThreeDSeismicModel Model;
    jif3D::TomographyData Data;
    Model.SetCellSize(grid.h, ncells, ncells, ncells);
    Data.AddMeasurementPoint(geo.x[1], geo.y[1], geo.z[1]);
    Data.AddSource(geo.x[0], geo.y[0], geo.z[0]);
    Data.AddMeasurementPoint(geo.x[0], geo.y[0], geo.z[0]);
    Data.AddSource(geo.x[1], geo.y[1], geo.z[1]);
    Data.AddMeasurementConfiguration(0, 0);
    Data.AddMeasurementConfiguration(1, 1);

    std::fill_n(Model.SetSlownesses().origin(), ncells * ncells * ncells, 0.1);

    jif3D::rvec time(Calculator.Calculate(Model, Data));
    jif3D::rvec Error(time.size());
    std::generate(Error.begin(),Error.end(),rand);
    //compare the C++ result to the C-result
    //the C++ class adds a low velocity layer at the bottom and the top
    //so the models are not identical, still the influence of these layers should be small
    BOOST_CHECK_CLOSE(data.tcalc[0], time(0),0.2);
    BOOST_CHECK_CLOSE(data.tcalc[1], time(1),0.2);
    //check the agreement between raypath and modelling for the calculator object
    totallength = 0.0;
    for (size_t i = 0; i < Calculator.GetRayPath()[0].nray + 1; ++i)
      {
        totallength += Calculator.GetRayPath()[0].len[i];
      }

    //check that the difference between calculated and exact solution
    //is smaller than the time it takes to go through one grid cell
    BOOST_CHECK(std::abs(dist-totallength * grid.h) < slow * grid.h);
    //now check that saving and restoring works
    Data.SetDataAndErrors(std::vector<double>(time.begin(),time.end()),std::vector<double>(Error.begin(),Error.end()));
    Data.WriteNetCDF("tt.nc");

    jif3D::TomographyData ReadData;
    ReadData.ReadNetCDF("tt.nc");
    BOOST_CHECK(std::equal(Data.GetData().begin(),Data.GetData().end(),ReadData.GetData().begin()));
    BOOST_CHECK(std::equal(Data.GetErrors().begin(),Data.GetErrors().end(),ReadData.GetErrors().begin()));
    BOOST_CHECK(std::equal(Data.GetSourcePosX().begin(),Data.GetSourcePosX().end(),ReadData.GetSourcePosX().begin()));
    BOOST_CHECK(std::equal(Data.GetSourcePosY().begin(),Data.GetSourcePosY().end(),ReadData.GetSourcePosY().begin()));
    BOOST_CHECK(std::equal(Data.GetSourcePosZ().begin(),Data.GetSourcePosZ().end(),ReadData.GetSourcePosZ().begin()));
    BOOST_CHECK(std::equal(Data.GetMeasPosX().begin(),Data.GetMeasPosX().end(),ReadData.GetMeasPosX().begin()));
    BOOST_CHECK(std::equal(Data.GetMeasPosY().begin(),Data.GetMeasPosY().end(),ReadData.GetMeasPosY().begin()));
    BOOST_CHECK(std::equal(Data.GetMeasPosZ().begin(),Data.GetMeasPosZ().end(),ReadData.GetMeasPosZ().begin()));
  }
BOOST_AUTO_TEST_SUITE_END()
