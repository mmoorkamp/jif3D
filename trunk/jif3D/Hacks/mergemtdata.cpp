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
#include "../ModelBase/VTKTools.h"
#include "../MT/ReadWriteImpedances.h"
#include "../MT/MTEquations.h"

int main()
  {
    std::string ncfilename1 = jif3D::AskFilename("Name of first netcdf file: ");
    std::string ncfilename2 = jif3D::AskFilename("Name of second netcdf file: ");

    jif3D::rvec Impedances1, Errors1;
    std::vector<double> Frequencies1, StatX1, StatY1, StatZ1, C1;
    jif3D::ReadImpedancesFromNetCDF(ncfilename1, Frequencies1, StatX1, StatY1, StatZ1,
        Impedances1, Errors1, C1);

    jif3D::rvec Impedances2, Errors2;
    std::vector<double> Frequencies2, StatX2, StatY2, StatZ2, C2;
    jif3D::ReadImpedancesFromNetCDF(ncfilename2, Frequencies2, StatX2, StatY2, StatZ2,
        Impedances2, Errors2, C2);

    jif3D::rvec Impedances(Impedances1.size() + Impedances2.size()), Errors(
        Errors1.size() + Errors2.size());
    std::vector<double> Frequencies(Frequencies1.size()), StatX(
        StatX1.size() + StatX2.size()), StatY(StatY1.size() + StatY2.size()), StatZ(
        StatZ1.size() + StatZ2.size()), C(C1.size() + C2.size());

    std::copy(StatX1.begin(), StatX1.end(), StatX.begin());
    std::copy(StatX2.begin(), StatX2.end(), StatX.begin() + StatX1.size());

    std::copy(StatY1.begin(), StatY1.end(), StatY.begin());
    std::copy(StatY2.begin(), StatY2.end(), StatY.begin() + StatY1.size());

    std::copy(StatZ1.begin(), StatZ1.end(), StatZ.begin());
    std::copy(StatZ2.begin(), StatZ2.end(), StatZ.begin() + StatZ1.size());

    std::copy(C1.begin(), C1.end(), C.begin());
    std::copy(C2.begin(), C2.end(), C.begin() + C1.size());


    for (size_t i = 0; i < Frequencies.size(); ++i)
      {
        size_t shift1 = StatX1.size() * 8;
        size_t shift2 = StatX2.size() * 8;
        std::copy(Impedances1.begin() + shift1 * i,
            Impedances1.begin() + shift1 * (i + 1),
            Impedances.begin() + (shift1 + shift2) * i);
        std::copy(Impedances2.begin() + shift2 * i,
            Impedances2.begin() + shift2 * (i + 1),
            Impedances.begin() + (shift1 + shift2) * i + shift1);

        std::copy(Errors1.begin() + shift1 * i,
            Errors1.begin() + shift1 * (i + 1),
            Errors.begin() + (shift1 + shift2) * i);
        std::copy(Errors2.begin() + shift2 * i,
            Errors2.begin() + shift2 * (i + 1),
            Errors.begin() + (shift1 + shift2) * i + shift1);
      }

    jif3D::WriteImpedancesToNetCDF(ncfilename1 + ".merged.nc",Frequencies1,StatX,StatY,StatZ,Impedances,Errors,C);
  }
