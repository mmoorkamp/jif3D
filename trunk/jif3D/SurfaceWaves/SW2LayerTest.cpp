/*
 * SW2LayerTest.cpp
 *
 *  Created on: 10.04.2020
 *      Author: bmw
 */

#include "../SurfaceWaves/SurfaceWaveCalculator.h"
#include <fstream>

int main()
  {
    jif3D::SurfaceWaveModel TstMod;

    std::vector<double> depth =
      { 2.0, 4.3, 6.8, 9.6, 12.8, 20.0 };
    //{ 35000.0, 75000.0, 100000.0 };

    std::vector<double> alpha =
     { 650.0, 750.0, 1400.0, 1800.0, 2150.0, 2800.0 };
     std::vector<double> beta =
     { 194.0, 270.0, 367.0, 485.0, 603.0, 740.0 };
     std::vector<double> rho =
     { 1820.0, 1860.0, 1910.0, 1960.0, 2020.0, 2090.0 };
    /*std::vector<double> alpha =
     { 6000.0, 8000.0 };
     std::vector<double> beta =
     { 3500.0, 4500.0 };
     std::vector<double> rho =
     { 2700.0, 3300.0 };
    std::vector<double> alpha =
      { 6000.0, 7000.0, 8000.0 };
    std::vector<double> beta =
      { 3500.0, 4000.0, 4500.0 };
    std::vector<double> rho =
      { 2700.0, 3000.0, 3300.0 };*/

    std::vector<double> PhaseTravelTimesErrors(1, 1.0);
    std::vector<double> Period =
    //{ 20.0, 30.0, 40.0 };
          { 1.0 / 30.0, 1.0 / 25.0, 1.0 / 20.0, 1.0 / 15.0, 1.0 / 10.0, 1.0 / 5.0 };
    /*std::vector<double> Period(80);
     for (int T = 0; T < 80; ++T)
     {
     Period[T] = (double) T + 1.0;
     }*/

    jif3D::SurfaceWaveCalculator TstCalc;
    TstCalc.set_data_err(PhaseTravelTimesErrors);

    std::ofstream outfile("surf.out");
    const size_t nfreqs = Period.size();
    std::vector<double> w(nfreqs);
    for (size_t i = 0; i < nfreqs; ++i)
      w[i] = 2 * M_PI / (Period[i]);

    for (size_t i = 0; i < nfreqs; ++i)
      {
        jif3D::SurfaceWaveCalculator::Surf1DResult Result = TstCalc.CalcSurf1D(w, i, rho,
            beta, alpha, depth);
        outfile << 2 * M_PI / w[i] << " " << Result.c << " " << Result.rc;
        for (int n = 0; n < depth.size(); ++n)
          {
            outfile << " " << Result.dcdvs[n] << " " << Result.dcdvp[n] << " "
                << Result.dcdrho[n];
          }
        outfile << std::endl;
      }

    return 0;
  }

