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

    std::vector<double> Impedances, Errors;
    std::vector<double> Frequencies, StatX, StatY, StatZ, C;
    std::vector<std::string> Names;

    jif3D::ReadImpedancesFromNetCDF(datafilename, Frequencies, StatX, StatY, StatZ,
        Impedances, Errors, C, Names);

    std::vector<double> SynthImpedances, SynthErrors;
    std::vector<double> SynthFrequencies, SynthStatX, SynthStatY, SynthStatZ, SynthC;
    std::vector<std::string> SynthNames;

    std::string synthfilename = jif3D::AskFilename("Name of synthetic file: ");
    jif3D::ReadImpedancesFromNetCDF(synthfilename, SynthFrequencies, SynthStatX,
        SynthStatY, SynthStatZ, SynthImpedances, SynthErrors, SynthC,SynthNames);

    size_t nshift = 0;
    for (size_t i = 0; i < Errors.size() - 1; i += 2)
      {
        if (Errors.at(i)
            > 1.0 * sqrt(jif3D::pow2(Impedances.at(i)) + jif3D::pow2(Impedances.at(i + 1))))
          {
            Impedances.at(i) = SynthImpedances.at(i);
            Impedances.at(i + 1) = SynthImpedances.at(i + 1);
            Errors.at(i) = 10.0 * SynthImpedances.at(i);
            Errors.at(i+1) = 10.0 * SynthImpedances.at(i+1);
            ++nshift;
          }

      }
    std::cout << "Modified: " << nshift * 2 << " out of " << Impedances.size() << " data "
        << std::endl;
    jif3D::WriteImpedancesToNetCDF(datafilename + ".shf.nc", Frequencies, StatX, StatY,
        StatZ, Impedances, Errors, C, Names);
  }
