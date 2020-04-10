/*
 * SW2LayerTest.cpp
 *
 *  Created on: 10.04.2020
 *      Author: bmw
 */

#include "../SurfaceWaves/SurfaceWaveCalculator.h"

int main()
  {
    jif3D::SurfaceWaveModel TstMod;

    std::vector<double> xcoords_m =
      { 0.0, 10000.0, 20000.0 };
    std::vector<double> ycoords_m =
      { 500000.0, 500000.0 + 10000.0, 500000.0 + 20000.0 };
    std::vector<double> zcoords_m =
      { 0.0, 35000.0, 100000.0 };

    std::vector<double> alpha =
      { 6000.0, 8000.0 };
    std::vector<double> beta =
      { 3500.0, 4500.0 };
    std::vector<double> rho =
      { 2700.0, 3300.0 };

    TstMod.SetCellCoords(xcoords_m, ycoords_m, zcoords_m);
    for (size_t i = 0; i < TstMod.GetData().num_elements(); ++i)
      {
        double Vs = beta[i % (zcoords_m.size() - 1)];
        double Vp = alpha[i % (zcoords_m.size() - 1)];
        double Density = rho[i % (zcoords_m.size() - 1)];

        TstMod.SetData().origin()[i] = Vs;
        TstMod.SetVp().origin()[i] = Vp;
        TstMod.SetDens().origin()[i] = Density;
      }

    jif3D::SurfaceWaveData TstDat;

    std::vector<double> xcoords_d =
      { 5000.0, 5000.0 };
    std::vector<double> ycoords_d =
      { 500000.0 + 5000.0, 500000.0 + 15000.0 };
    std::vector<double> zcoords_d =
      { 0.0, 0.0 };

    std::vector<double> PhaseTravelTimes(1, 1.0);
    std::vector<double> PhaseTravelTimesErrors(1, 1.0);
    std::vector<double> Period =
      { 20.0 };

    std::vector<double> EventPosX =
      { 0.0 };
    std::vector<double> EventPosY =
      { 0.0 };
    std::vector<double> EventPosZ =
      { 0.0 };
    std::vector<double> EventStationCombination =
      { 0 };
    std::vector<double> StationCombination =
      { 0, 1 };

    TstDat.SetMeasurementPoints(xcoords_d, ycoords_d, zcoords_d);
    TstDat.SetDataAndErrors(PhaseTravelTimes, PhaseTravelTimesErrors);
    TstDat.SetEventPositions(EventPosX, EventPosY, EventPosZ);
    TstDat.SetPeriods(Period);
    TstDat.SetEventStatCmb(EventStationCombination);
    TstDat.SetStatComb(StationCombination);
    TstDat.SetDummy(-999.9);
    TstDat.SetLonCentr(-123.0);

    jif3D::SurfaceWaveCalculator TstCalc;
    TstCalc.set_data_err(PhaseTravelTimesErrors);
    jif3D::rvec ModTravelTime(TstCalc.Calculate(TstMod, TstDat));

    return 0;
  }

