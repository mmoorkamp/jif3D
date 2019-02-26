//============================================================================
// Name        : qqfit.cpp
// Author      : Dec 6, 2013
// Version     :
// Copyright   : 2013, mmoorkamp
//============================================================================

#include <iostream>
#include <fstream>
#include <string>
#include <boost/program_options.hpp>
#include "../Global/FileUtil.h"
#include "../Global/NumUtil.h"
#include "../Global/FatalException.h"
#include "../MT/ReadWriteImpedances.h"
#include "../MT/MTEquations.h"

int main()
  {
    std::string datafilename = jif3D::AskFilename("Name of data file: ");

    jif3D::rvec Impedances, Errors;
    std::vector<double> Frequencies, StatX, StatY, StatZ, C;
    jif3D::ReadImpedancesFromNetCDF(datafilename, Frequencies, StatX, StatY, StatZ,
        Impedances, Errors, C);

    jif3D::rvec SynthImpedances, SynthErrors;
    std::vector<double> SynthFrequencies, SynthStatX, SynthStatY, SynthStatZ, SynthC;
    std::string synthfilename = jif3D::AskFilename("Name of synthetic file: ");
    jif3D::ReadImpedancesFromNetCDF(synthfilename, SynthFrequencies, SynthStatX,
        SynthStatY, SynthStatZ, SynthImpedances, SynthErrors, SynthC);

    size_t nshift = 0;
    for (size_t i = 0; i < Errors.size() - 1; i += 2)
      {
        if (Errors(i)
            > 0.5 * sqrt(jif3D::pow2(Impedances(i)) + jif3D::pow2(Impedances(i + 1))))
          {
            Impedances(i) = SynthImpedances(i);
            Impedances(i + 1) = SynthImpedances(i + 1);
            ++nshift;
          }

      }
    std::cout << "Modified: " << nshift * 2 << " out of " << Impedances.size() << " data "
        << std::endl;
    jif3D::WriteImpedancesToNetCDF(datafilename + ".shf.nc", Frequencies, StatX, StatY,
        StatZ, Impedances, Errors, C);
  }
