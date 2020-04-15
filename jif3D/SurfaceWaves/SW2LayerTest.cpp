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
      { 35000.0, 100000.0 };

    std::vector<double> alpha =
      { 6000.0, 8000.0 };
    std::vector<double> beta =
      { 3500.0, 4500.0 };
    std::vector<double> rho =
      { 2700.0, 3300.0 };

    std::vector<double> PhaseTravelTimesErrors(1, 1.0);
    std::vector<double> Period =
      { 20.0 };

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
        outfile << 2 * M_PI / w[i] << " " << Result.c << " " << Result.rc << " "
            << Result.dcdvp[0] << " " << Result.dcdvp[1] << std::endl;
      }

    return 0;
  }

