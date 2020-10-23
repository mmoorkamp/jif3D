/*
 * test_surfacewave_gradient.cpp
 *
 *  Created on: 20 Nov 2019
 *      Author: bweise
 */

#define BOOST_TEST_MODULE SurfaceWaveDisp test
#define BOOST_TEST_MAIN ...
#include "../Global/Jif3DTesting.h"
#include <boost/test/included/unit_test.hpp>
#include <stdlib.h>
#include <time.h>
#include <vector>
#include "../SurfaceWaves/SurfaceWaveCalculator.h"
#include "../SurfaceWaves/SurfaceWaveFunctions.h"
#include "../Inversion/ModelTransforms.h"
#include "../Inversion/ThreeDModelObjective.h"
#include <boost/test/tools/floating_point_comparison.hpp>
#include <boost/random/linear_congruential.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>
#include <GeographicLib/TransverseMercatorExact.hpp>
#include <GeographicLib/Constants.hpp>

BOOST_AUTO_TEST_SUITE (SW_Dispersion_Test_Suite)

    BOOST_AUTO_TEST_CASE (oned_test)
      {
        jif3D::SurfaceWaveCalculator Calculator;
        std::vector<double> alpha =
          { 6000.0, 8000.0, 8000.0 };
        std::vector<double> beta =
          { 3500.0, 4500.0, 4500.0 };
        std::vector<double> rho =
          { 2700.0, 3300.0, 3300.0 };
        //std::vector<double> T =
        //  { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 120, 140,
        //      160, 180, 200 };
        std::vector<double> T =
          { 1 };
        std::vector<double> w = jif3D::T2w(T);

        std::vector<double> zcoords_m =
        //{ 0.0, 2.0, 4.3, 6.8, 9.6, 12.8, 20.0 };
              { 0.0, 5000.0, 10000.0, 15000.0 };
        const std::vector<double> depth(zcoords_m.begin() + 1, zcoords_m.end());
        std::vector<double> oneddisp(T.size());

        std::vector<double> statlat =
          { 45.0, 46.0 };
        std::vector<double> statlon =
          { -123.0, -122.0 };
        GeographicLib::TransverseMercatorExact proj(GeographicLib::Constants::WGS84_a(),
            GeographicLib::Constants::WGS84_f(), GeographicLib::Constants::UTM_k0());

        std::vector<double> ycoords(statlat.size()), xcoords(statlat.size()), zcoords(
            statlat.size());
        double lon_centr = -123.0;
        double false_east = 500000.0;
        for (int i = 0; i < statlat.size(); ++i)
          {
            double x, y;
            proj.Forward(lon_centr, statlat.at(i), statlon.at(i), y, x);
            xcoords.at(i) = x;
            ycoords.at(i) = y + false_east;
            zcoords.at(i) = 0.0;
            std::cout << "Stat X: " << xcoords.at(i) << " Stat Y:" << ycoords.at(i)
                << std::endl;
          }

        const int start = -2;
        const int ncells = 20;
        const double step = 10000.0;
        std::vector<double> xcoords_m(ncells), ycoords_m(ncells);
        for (int i = start; i < ncells + start; ++i)
          {
            xcoords_m.at(i - start) = xcoords.at(0) + i * step;
            ycoords_m.at(i - start) = ycoords.at(0) + i * step;
          }

        jif3D::SurfaceWaveModel TomoModel;
        const int zsize = zcoords_m.size();

        TomoModel.SetCellCoords(xcoords_m, ycoords_m, zcoords_m);

        auto exten =
            boost::extents[xcoords_m.size() - 1][ycoords_m.size() - 1][zcoords_m.size()
                - 1];
        jif3D::SurfaceWaveModel::t3DModelData Velocity(exten), Vp(exten), Density(exten);
        for (size_t i = 0; i < TomoModel.GetData().num_elements(); ++i)
          {

            Velocity.origin()[i] = beta[i % (zsize - 1)];
            Vp.origin()[i] = alpha[i % (zsize - 1)];
            Density.origin()[i] = rho[i % (zsize - 1)];

            TomoModel.SetData() = Velocity;
            TomoModel.SetVp(Vp);
            TomoModel.SetDens(Density);
          }

        jif3D::SurfaceWaveData SWData;
        SWData.SetMeasurementPoints(xcoords, ycoords, zcoords);

        SWData.SetPeriods(T);
        SWData.SetLonCentr(lon_centr);
        std::vector<int> stat1 =
          { 0 };
        std::vector<int> stat2 =
          { 1 };
        SWData.SetStationPairs(stat1, stat2);
        std::vector<int> NDataPerT(T.size(), 1);
        SWData.SetNDataPerT(NDataPerT);
        std::vector<int> PairIndex(T.size(), 0);

        std::vector<int> EventIndex(T.size(), 0);

        std::vector<int> PeriodIndex(T.size());
        std::iota(PeriodIndex.begin(), PeriodIndex.end(), 0);

        std::vector<double> dtp(PairIndex.size(), 1.0);
        std::vector<double> err(dtp.size(), 1.0);
        SWData.SetIndexMap(PairIndex, EventIndex, PeriodIndex);
        SWData.SetDataAndErrors(dtp, err);

//Calculator.set_data_err(err);
        double dist = 200000.0;
        double xcent = (xcoords.at(0) + xcoords.at(1)) / 2.0;
        double ycent = (ycoords.at(0) + ycoords.at(1)) / 2.0;

        std::ofstream anglefile("disp_angle.out");
        for (double angle = 0; angle < 360; angle += 5.0)
          {
            double eventx = xcent + dist * std::cos(angle / 180.0 * M_PI);
            double eventy = ycent + dist * std::sin(angle / 180.0 * M_PI);
            double elat, elon;
            proj.Reverse(lon_centr, eventy - false_east, eventx, elat, elon);

            std::vector<double> eventlat =
              { elat };
            std::vector<double> eventlon =
              { elon };
            std::vector<double> epz =
              { 100000.0 };
            SWData.SetEventPositions(eventlat, eventlon, epz);
            jif3D::rvec ObservedTimes(Calculator.Calculate(TomoModel, SWData));
            anglefile << angle << " " << ObservedTimes(0) << std::endl;
          }

        jif3D::rvec ObservedTimes(Calculator.Calculate(TomoModel, SWData));

        double statdist = std::sqrt(
            std::pow(xcoords.at(0) - xcoords.at(1), 2)
                + std::pow(ycoords.at(0) - ycoords.at(1), 2));
        std::cout << "Dist: " << dist << std::endl;
        std::ofstream dispfile("disp.out");
        for (int index = 0; index < T.size(); ++index)
          {
            jif3D::SurfaceWaveCalculator::Surf1DResult Result = Calculator.CalcSurf1D(w,
                index, rho, beta, alpha, depth);
            oneddisp.at(index) = Result.c;
            dispfile << T.at(index) << " " << Result.c << " " << ObservedTimes(index)
                << " " << statdist / ObservedTimes(index) << std::endl;
          }

      }

    BOOST_AUTO_TEST_SUITE_END()
