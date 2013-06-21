//============================================================================
// Name        : writemtt.cpp
// Author      : Feb 21, 2011
// Version     :
// Copyright   : 2011, mmoorkamp
//============================================================================

#include <iostream>
#include <fstream>
#include <string>
#include <boost/program_options.hpp>
#include "../Global/FileUtil.h"
#include "../Global/FatalException.h"
#include "../MT/ReadWriteImpedances.h"
#include "../MT/MTEquations.h"

int main()
  {
    std::string ncfilename = jif3D::AskFilename("Name of netcdf file: ");

    jif3D::rvec Impedances, Errors;
    std::vector<double> Frequencies, StatX, StatY, StatZ;
    jif3D::ReadImpedancesFromNetCDF(ncfilename, Frequencies, StatX, StatY,
        StatZ, Impedances, Errors);
    jif3D::WriteImpedancesToMtt(ncfilename, Frequencies, Impedances, Errors);

    std::ofstream Zxx_rho("Zxx_rho.xyz");
    std::ofstream Zxy_rho("Zxy_rho.xyz");
    std::ofstream Zyx_rho("Zyx_rho.xyz");
    std::ofstream Zyy_rho("Zyy_rho.xyz");

    std::ofstream Zxx_phi("Zxx_phi.xyz");
    std::ofstream Zxy_phi("Zxy_phi.xyz");
    std::ofstream Zyx_phi("Zyx_phi.xyz");
    std::ofstream Zyy_phi("Zyy_phi.xyz");
    std::ofstream Posfile("station.pos");

    const size_t nimp = Impedances.size();
    const size_t nval = nimp / 8;
    const size_t nstats = StatX.size();
    for (size_t i = 0; i < nstats; ++i)
      {
        Posfile << StatX.at(i) << " " << StatY.at(i) << " " << StatZ.at(i) << std::endl;
      }
    for (size_t i = 0; i < nval; ++i)
      {
        size_t freqindex = i / nstats;
        size_t statindex = i % nstats;
        Zxx_rho << statindex << " " << Frequencies.at(freqindex) << " "
            << jif3D::AppRes(std::complex<double>(Impedances(i * 8), Impedances(
                i * 8 + 1)), Frequencies.at(freqindex)) << "\n";

        Zxy_rho << statindex << " " << Frequencies.at(freqindex) << " "
            << jif3D::AppRes(std::complex<double>(Impedances(i * 8 + 2),
                Impedances(i * 8 + 3)), Frequencies.at(freqindex)) << "\n";

        Zyx_rho << statindex << " " << Frequencies.at(freqindex) << " "
            << jif3D::AppRes(std::complex<double>(Impedances(i * 8 + 4),
                Impedances(i * 8 + 5)), Frequencies.at(freqindex)) << "\n";

        Zyy_rho << statindex << " " << Frequencies.at(freqindex) << " "
            << jif3D::AppRes(std::complex<double>(Impedances(i * 8 + 6),
                Impedances(i * 8 + 7)), Frequencies.at(freqindex)) << "\n";

        Zxx_phi << statindex << " " << Frequencies.at(freqindex) << " "
            << jif3D::ImpedancePhase(std::complex<double>(Impedances(i * 8),
                Impedances(i * 8 + 1))) << "\n";

        Zxy_phi << statindex << " " << Frequencies.at(freqindex) << " "
            << jif3D::ImpedancePhase(std::complex<double>(Impedances(i * 8 + 2),
                Impedances(i * 8 + 3))) << "\n";

        Zyx_phi << statindex << " " << Frequencies.at(freqindex) << " "
            << jif3D::ImpedancePhase(std::complex<double>(Impedances(i * 8 + 4),
                Impedances(i * 8 + 5))) << "\n";

        Zyy_phi << statindex << " " << Frequencies.at(freqindex) << " "
            << jif3D::ImpedancePhase(std::complex<double>(Impedances(i * 8 + 6),
                Impedances(i * 8 + 7))) << "\n";
      }

  }
