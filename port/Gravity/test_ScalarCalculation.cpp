//============================================================================
// Name        : test_ScalarCalculation.cpp
// Author      : Feb 11, 2009
// Version     :
// Copyright   : 2009, mmoorkamp
//============================================================================

#define BOOST_TEST_MODULE ThreeDGravityModel test
#define BOOST_TEST_MAIN ...
#include <boost/test/included/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/assign/std/vector.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/numeric/conversion/cast.hpp>
#include <boost/math/constants/constants.hpp>

#include "test_common.h"
#include "BasicGravElements.h"
#include "../GravMag/MinMemGravMagCalculator.h"
#include "../GravMag/FullSensitivityGravMagCalculator.h"
#include "ScalarOMPGravityImp.h"
#include "ThreeDGravityFactory.h"
#include "../Global/Jif3DPlatformHelper.h"

using namespace boost::assign;

BOOST_AUTO_TEST_SUITE( ScalarGravity_Test_Suite )

//compare calculations for a single prims with known values
    BOOST_AUTO_TEST_CASE (box_gravity_calc_test)
      {
        //the gravity in the center of the cube should be 0
        //the density is 1 in all calculations here so do not have to multiply the term
        BOOST_CHECK(
            fabs(
                jif3D::CalcGravBoxTerm(0.0, 0.0, 0.0, -10.0, -10.0, -10.0, 20.0, 20.0,
                    20.0)) < std::numeric_limits<double>::epsilon());
        //compare with the reported value of li and chouteau within a precision of 0.1%
        double topofcube = jif3D::CalcGravBoxTerm(0.0, 0.0, -10.0, -10.0, -10.0, -10.0,
            20.0, 20.0, 20.0);
        BOOST_CHECK_CLOSE(topofcube, 3.46426e-6, 0.1);
        double bottomofcube = jif3D::CalcGravBoxTerm(0.0, 0.0, 10.0, -10.0, -10.0, -10.0,
            20.0, 20.0, 20.0);
        //check symmetry of results
        BOOST_CHECK_CLOSE(topofcube, -bottomofcube,
            std::numeric_limits<float>::epsilon());
        //check translation invariance
        double shifted = jif3D::CalcGravBoxTerm(10.0, 10.0, 0.0, 0.0, 0.0, 0.0, 20.0,
            20.0, 20.0);
        BOOST_CHECK_CLOSE(topofcube, shifted, std::numeric_limits<float>::epsilon());

        double away1 = jif3D::CalcGravBoxTerm(1e5, 10.0, 10.0, -10.0, -10.0, -10.0, 20.0,
            20.0, 20.0);
        BOOST_CHECK(away1 < std::numeric_limits<double>::epsilon());
        double away2 = jif3D::CalcGravBoxTerm(100.0, 100.0, 10.0, -10.0, -10.0, -10.0,
            20.0, 20.0, 20.0);
        BOOST_CHECK_CLOSE(away2, -1.873178e-9, 0.1);
        //check also for a very large cube
        BOOST_CHECK(
            fabs(jif3D::CalcGravBoxTerm(0.0, 0.0, 0.0, -1e6, -1e6, -1e6, 2e6, 2e6, 2e6))
                < std::numeric_limits<double>::epsilon());
      }

    //check that for a box the results are independent of the discretization
    BOOST_AUTO_TEST_CASE(model_gravity_boxcomp_test)
      {
        const double measx = 9.0;
        const double measy = 8.0;
        const double measz = -0.1;
        //again density is 1
        double boxtopofcube = jif3D::CalcGravBoxTerm(measx, measy, measz, 0.0, 0.0, 0.0,
            20.0, 20.0, 20.0);
        jif3D::GravimetryMatrix tensorbox = jif3D::CalcTensorBoxTerm(measx, measy, measz,
            0.0, 0.0, 0.0, 20.0, 20.0, 20.0);
        jif3D::ThreeDGravityModel GravityTest;
        //create a model of 10x10x10 cells with 2m length in each dimension
        const size_t ncells = 10;
        const double cellsize = 2.0;
        GravityTest.SetMeshSize(ncells,ncells,ncells);
        for (size_t i = 0; i < ncells; ++i)
          {
            GravityTest.SetXCellSizes()[i] = cellsize;
            GravityTest.SetYCellSizes()[i] = cellsize;
            GravityTest.SetZCellSizes()[i] = cellsize;
          }
        jif3D::rvec DensityVector(ncells * ncells * ncells, 1.0);
        std::fill_n(GravityTest.SetDensities().origin(), ncells * ncells * ncells, 1.0);

        GravityTest.AddMeasurementPoint(measx, measy, measz);
        typedef typename jif3D::FullSensitivityGravMagCalculator<jif3D::ThreeDGravityModel> CalculatorType;
        boost::shared_ptr<CalculatorType> TensorCalculator(
            jif3D::CreateGravityCalculator<CalculatorType>::MakeTensor());
        boost::shared_ptr<CalculatorType> ScalarCalculator(
            jif3D::CreateGravityCalculator<CalculatorType>::MakeScalar());
        jif3D::rvec gravmeas(ScalarCalculator->Calculate(GravityTest));
        jif3D::rvec tensmeas(TensorCalculator->Calculate(GravityTest));
        double gridcube = gravmeas[0];

        GravityTest.WriteNetCDF("cube.nc");
        // Check that FTG calculation does not depend on discretization
        for (size_t i = 0; i < 9; ++i)
          {
            BOOST_CHECK_CLOSE(tensmeas(i), tensorbox.data()[i],
                std::numeric_limits<float>::epsilon());
          }

        //check some tensor properties
        BOOST_CHECK_CLOSE(tensmeas(0) + tensmeas(4), -tensmeas(8),
            std::numeric_limits<float>::epsilon());
        BOOST_CHECK_CLOSE(boxtopofcube, gridcube, std::numeric_limits<float>::epsilon());
        //Check scalar sensitivity calculation
        jif3D::rmat ScalarSensitivities(ScalarCalculator->GetSensitivities());
        jif3D::rvec SensValues(
            boost::numeric::ublas::prec_prod(ScalarSensitivities, DensityVector));
        BOOST_CHECK_CLOSE(gridcube, SensValues(0), std::numeric_limits<float>::epsilon());
      }

    //check that the 1D background calculation gives correct results
    BOOST_AUTO_TEST_CASE(background_test)
      {
        jif3D::ThreeDGravityModel GravityTest;
        double analytic = 2 * boost::math::constants::pi<double>() * jif3D::Grav_const
            * (500.0 + 5.0 * 4500.0);
        std::vector<double> bg_dens, bg_thick;
        const size_t nmeas = 10;
        for (size_t i = 0; i < nmeas; ++i)
          GravityTest.AddMeasurementPoint(rand() % 50000 + 2e4, rand() % 50000 + 2e4,
              0.0);
        bg_dens += 1.0, 1.0, 5.0, 5.0;
        bg_thick += 200.0, 300.0, 3500.0, 1000.0;
        GravityTest.SetBackgroundDensities(bg_dens);
        GravityTest.SetBackgroundThicknesses(bg_thick);
        typedef typename jif3D::FullSensitivityGravMagCalculator<jif3D::ThreeDGravityModel> CalculatorType;
        boost::shared_ptr<CalculatorType> ScalarCalculator(
            jif3D::CreateGravityCalculator<CalculatorType>::MakeScalar());
        jif3D::rvec gravmeas(ScalarCalculator->Calculate(GravityTest));
        for (size_t i = 0; i < nmeas; ++i)
          {
            BOOST_CHECK_CLOSE(analytic, gravmeas[i], 0.01);
            BOOST_CHECK_CLOSE(gravmeas[0], gravmeas[i], 0.01);
          }
      }

    //compare the 3D solution with the analytic solution for a 1D layered structure
    BOOST_AUTO_TEST_CASE(model_gravity_1danaly_test)
      {
        jif3D::ThreeDGravityModel GravityTest;
        //create a model with a random number of cells in each direction
        const size_t nhorcells = rand() % 50 + 10;
        const size_t nzcells = 10;
        GravityTest.SetMeshSize(nhorcells, nhorcells, nzcells);
        jif3D::rvec DensityVector(nhorcells * nhorcells * nzcells + 4); // 4 background layers
        for (size_t i = 0; i < nhorcells; ++i) // set the values of the inner cells

          {
            GravityTest.SetXCellSizes()[i] = rand() % 10000 + 1000;
            GravityTest.SetYCellSizes()[i] = rand() % 10000 + 1000;
          }
        for (size_t i = 0; i < nzcells; ++i)
          {
            GravityTest.SetZCellSizes()[i] = 500;
          }


        for (size_t i = 0; i < nhorcells; ++i)
          for (size_t j = 0; j < nhorcells; ++j)
            {
              GravityTest.SetDensities()[i][j][0] = 1.0;
              DensityVector(i * nhorcells * nzcells + j * nzcells) = 1.0;
            }
        for (size_t i = 0; i < nhorcells; ++i)
          for (size_t j = 0; j < nhorcells; ++j)
            for (size_t k = 1; k < nzcells; ++k)
              {
                GravityTest.SetDensities()[i][j][k] = 5.0;
                DensityVector(i * nhorcells * nzcells + j * nzcells + k) = 5.0;
              }

        const size_t nmeas = 10;
        for (size_t i = 0; i < nmeas; ++i)
          GravityTest.AddMeasurementPoint(rand() % 50000 + 2e4, rand() % 50000 + 2e4,
              -1.0);
        std::vector<double> bg_dens, bg_thick;
        bg_dens += 1.0, 1.0, 5.0, 5.0;
        bg_thick += 200.0, 300.0, 3500.0, 2000.0;

        DensityVector(nhorcells * nhorcells * nzcells) = 1.0;
        DensityVector(nhorcells * nhorcells * nzcells + 1) = 1.0;
        DensityVector(nhorcells * nhorcells * nzcells + 2) = 5.0;
        DensityVector(nhorcells * nhorcells * nzcells + 3) = 5.0;
        GravityTest.SetBackgroundDensities(bg_dens);
        GravityTest.SetBackgroundThicknesses(bg_thick);

        typedef typename jif3D::FullSensitivityGravMagCalculator<jif3D::ThreeDGravityModel> CalculatorType;
        boost::shared_ptr<CalculatorType> TensorCalculator(
            jif3D::CreateGravityCalculator<CalculatorType>::MakeTensor());
        boost::shared_ptr<CalculatorType> ScalarCalculator(
            jif3D::CreateGravityCalculator<CalculatorType>::MakeScalar());
#ifdef HAVEGPU
        boost::shared_ptr<CalculatorType> CudaCalculator(jif3D::CreateGravityCalculator<CalculatorType>::MakeScalar(true));
        std::cout << "Calculating CUDA " << std::endl;
        jif3D::rvec cudameas(CudaCalculator->Calculate(GravityTest));
        std::cout << "Finished CUDA " << std::endl;
#endif
        jif3D::rvec scalarmeas(ScalarCalculator->Calculate(GravityTest));
        jif3D::rvec tensormeas(TensorCalculator->Calculate(GravityTest));

        GravityTest.WriteNetCDF("layer.nc");
        double analytic = 2 * boost::math::constants::pi<double>() * jif3D::Grav_const
            * (500.0 + 5.0 * 5500.0);
        jif3D::rmat ScalarSensitivities(ScalarCalculator->GetSensitivities());
        jif3D::rvec SensValues(prec_prod(ScalarSensitivities, DensityVector));
        //The second time the measurements will also be calculated from the sensitivities internally
        jif3D::rvec scalarmeas2(ScalarCalculator->Calculate(GravityTest));
        for (size_t i = 0; i < nmeas; ++i)
          {
            BOOST_CHECK_CLOSE(analytic, scalarmeas[i], 0.01);
            BOOST_CHECK_CLOSE(analytic, scalarmeas2[i], 0.001);
            BOOST_CHECK_CLOSE(analytic, SensValues[i], 0.01);
#ifdef HAVEGPU
            BOOST_CHECK_CLOSE(analytic, cudameas[i], 0.01);
#endif
            // The  tensorial elements should all be zero
            for (size_t j = 0; j < 9; ++j)
              {
                BOOST_CHECK(
                    fabs(tensormeas(i * 9 + j)) < std::numeric_limits<double>::epsilon());
              }
          }
      }

    //check whether calculation by the sensitivity matrix yields the same results
    //for scalar calculations
    BOOST_AUTO_TEST_CASE(scalar_caching_test)
      {
        jif3D::ThreeDGravityModel GravityTest;
        const size_t ncells = 10;
        const size_t nmeas = 10;
        MakeRandomModel(GravityTest, ncells, nmeas);
        typedef typename jif3D::FullSensitivityGravMagCalculator<jif3D::ThreeDGravityModel> CalculatorType;
        boost::shared_ptr<CalculatorType> ScalarCalculator(
            jif3D::CreateGravityCalculator<CalculatorType>::MakeScalar());

        //Calculate twice, once with normal calculation, once cached
        //and record the times
        boost::posix_time::ptime startfirst =
            boost::posix_time::microsec_clock::local_time();
        jif3D::rvec scalarmeas1(ScalarCalculator->Calculate(GravityTest));
        boost::posix_time::ptime endfirst =
            boost::posix_time::microsec_clock::local_time();
        jif3D::rvec scalarmeas2(ScalarCalculator->Calculate(GravityTest));
        boost::posix_time::ptime endsecond =
            boost::posix_time::microsec_clock::local_time();
        //the second time it should be much faster
        BOOST_CHECK(endfirst - startfirst > (endsecond - endfirst));
        for (size_t i = 0; i < nmeas; ++i)
          {
            BOOST_CHECK_CLOSE(scalarmeas1[i], scalarmeas2[i],
                std::numeric_limits<float>::epsilon());
          }
      }

    BOOST_AUTO_TEST_CASE (lqderivative_test)
      {
        jif3D::ThreeDGravityModel GravityTest;
        const size_t ncells = 10;
        const size_t nmeas = 10;
        MakeRandomModel(GravityTest, ncells, nmeas, false);
        typedef typename jif3D::FullSensitivityGravMagCalculator<jif3D::ThreeDGravityModel> CalculatorType;
        boost::shared_ptr<CalculatorType> ScalarCalculator(
            jif3D::CreateGravityCalculator<CalculatorType>::MakeScalar());
        jif3D::rvec Misfit(nmeas);
        std::generate(Misfit.begin(), Misfit.end(), jif3D::platform::drand48);
        jif3D::rvec Deriv(ScalarCalculator->LQDerivative(GravityTest, Misfit));
        ScalarCalculator->Calculate(GravityTest);
        jif3D::rvec Compare(
            2.0
                * boost::numeric::ublas::prec_prod(
                    ublas::trans(ScalarCalculator->GetSensitivities()), Misfit));
        //and test the caching, too
        jif3D::rvec Deriv2(ScalarCalculator->LQDerivative(GravityTest, Misfit));
        const size_t ngrid = GravityTest.GetDensities().num_elements();
        BOOST_CHECK(Deriv.size() == ngrid);
        BOOST_CHECK(Compare.size() == ngrid);
        BOOST_CHECK(Deriv2.size() == ngrid);
        for (size_t i = 0; i < ngrid; ++i)
          {
            BOOST_CHECK_CLOSE(Deriv(i), Compare(i),
                std::numeric_limits<float>::epsilon());
            BOOST_CHECK_CLOSE(Deriv2(i), Compare(i),
                std::numeric_limits<float>::epsilon());
          }
      }

#ifdef HAVEUBCCODE
    //compare with the ubc gravity forward code
    BOOST_AUTO_TEST_CASE(ubc_test)
      {
        //create a random model
        jif3D::ThreeDGravityModel GravityTest;
        const size_t nmeas = 10;
        const size_t ncells = 16;
        MakeRandomModel(GravityTest,ncells, nmeas);
        //and compute the results for our code
        typedef typename jif3D::MinMemGravMagCalculator<jif3D::ThreeDGravityModel> CalculatorType;
        boost::shared_ptr<CalculatorType> ScalarCalculator(jif3D::CreateGravityCalculator<CalculatorType>::MakeScalar());
        jif3D::rvec scalarmeas(
            ScalarCalculator->Calculate(GravityTest));
        //set the mesh for the ubc code, they have the x-axis in east direction
        std::ofstream meshfile("testmesh");
        meshfile << GravityTest.GetYCellSizes().size() << " ";
        meshfile << GravityTest.GetXCellSizes().size() << " ";
        meshfile << GravityTest.GetZCellSizes().size() << "\n";
        meshfile << " 0 0 0 \n";
        //so we have to swap x and y in the grid coordinate system
        std::copy(GravityTest.GetYCellSizes().begin(),
            GravityTest.GetYCellSizes().end(), std::ostream_iterator<double>(
                meshfile, " "));
        meshfile << "\n";
        std::copy(GravityTest.GetXCellSizes().begin(),
            GravityTest.GetXCellSizes().end(), std::ostream_iterator<double>(
                meshfile, " "));
        meshfile << "\n";
        std::copy(GravityTest.GetZCellSizes().begin(),
            GravityTest.GetZCellSizes().end(), std::ostream_iterator<double>(
                meshfile, " "));
        meshfile << "\n" << std::flush;
        //write the observation points
        //again we have to switch x and y and z is positive up
        std::ofstream obsfile("testobs");
        obsfile << nmeas << "\n";
        for (size_t i = 0; i < nmeas; ++i)
          {
            obsfile << GravityTest.GetMeasPosY().at(i) << " "
            << GravityTest.GetMeasPosX().at(i) << " "
            << -GravityTest.GetMeasPosZ().at(i) << "\n";
          }
        obsfile << std::flush;
        //write out the densities, we can write them in the same order we store them
        std::ofstream densfile("testdens");
        std::copy(GravityTest.GetDensities().origin(),
            GravityTest.GetDensities().origin()
            + GravityTest.GetDensities().num_elements(),
            std::ostream_iterator<double>(densfile, "\n"));
        densfile << std::flush;
        //call the ubc code using wine
        system("wine gzfor3d.exe testmesh testobs testdens");
        //read in the results
        std::ifstream resultfile("gzfor3d.grv");
        //the first line contains some extra information we don't want
        char dummy[255];
        resultfile.getline(dummy, 255);
        //read in all numbers in the file into a vector
        //this includes coordinate information etc.
        std::vector<double> ubcresults;
        std::copy(std::istream_iterator<double>(resultfile),
            std::istream_iterator<double>(), std::back_inserter(ubcresults));
        for (size_t i = 0; i < nmeas; ++i)
          {
            //compare our results with the right component of the vector
            //the tolerance of 0.1% is OK considering truncation effect from writing to ascii
            BOOST_CHECK_CLOSE(scalarmeas[i] * 100000.0,
                ubcresults[(i + 1) * 4 - 1], 0.1);
          }
      }
#endif

    BOOST_AUTO_TEST_SUITE_END()
