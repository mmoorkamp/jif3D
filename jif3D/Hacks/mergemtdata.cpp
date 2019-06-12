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
#include "../MT/MTData.h"

int main()
  {
    std::string ncfilename1 = jif3D::AskFilename("Name of first netcdf file: ");
    std::string ncfilename2 = jif3D::AskFilename("Name of second netcdf file: ");

    jif3D::MTData Data1, Data2;
    Data1.ReadNetCDF(ncfilename1);
    Data2.ReadNetCDF(ncfilename2);

    std::vector<double> Impedances(Data1.GetData().size() + Data1.GetData().size()), Errors(
        Data1.GetErrors().size() + Data1.GetErrors().size());
    std::vector<double> Frequencies(Data1.GetFrequencies().size()), StatX(
        Data1.GetMeasPosX().size() + Data2.GetMeasPosX().size()), StatY(
        Data1.GetMeasPosY().size() + Data2.GetMeasPosY().size()), StatZ(
        Data1.GetMeasPosZ().size() + Data2.GetMeasPosZ().size()), C(
        Data1.GetDistortion().size() + Data2.GetDistortion().size());

    std::copy(Data1.GetMeasPosX().begin(), Data1.GetMeasPosX().end(), StatX.begin());
    std::copy(Data2.GetMeasPosX().begin(), Data2.GetMeasPosX().end(),
        StatX.begin() + Data1.GetMeasPosX().size());

    std::copy(Data1.GetMeasPosY().begin(), Data1.GetMeasPosY().end(), StatY.begin());
    std::copy(Data2.GetMeasPosY().begin(), Data2.GetMeasPosY().end(),
        StatY.begin() + Data1.GetMeasPosY().size());

    std::copy(Data1.GetMeasPosZ().begin(), Data1.GetMeasPosZ().end(), StatZ.begin());
    std::copy(Data2.GetMeasPosZ().begin(), Data2.GetMeasPosZ().end(),
        StatZ.begin() + Data1.GetMeasPosZ().size());

    std::copy(Data1.GetDistortion().begin(), Data1.GetDistortion().end(), C.begin());
    std::copy(Data2.GetDistortion().begin(), Data2.GetDistortion().end(),
        C.begin() + Data1.GetDistortion().size());

    for (size_t i = 0; i < Frequencies.size(); ++i)
      {
        size_t shift1 = Data1.GetMeasPosX().size() * 8;
        size_t shift2 = Data2.GetMeasPosX().size() * 8;
        std::copy(Data1.GetData().begin() + shift1 * i,
            Data1.GetData().begin() + shift1 * (i + 1),
            Impedances.begin() + (shift1 + shift2) * i);
        std::copy(Data2.GetData().begin() + shift2 * i,
            Data2.GetData().begin() + shift2 * (i + 1),
            Impedances.begin() + (shift1 + shift2) * i + shift1);

        std::copy(Data1.GetErrors().begin() + shift1 * i, Data1.GetErrors().begin() + shift1 * (i + 1),
            Errors.begin() + (shift1 + shift2) * i);
        std::copy(Data2.GetErrors().begin() + shift2 * i, Data2.GetErrors().begin() + shift2 * (i + 1),
            Errors.begin() + (shift1 + shift2) * i + shift1);
      }

    jif3D::WriteImpedancesToNetCDF(ncfilename1 + ".merged.nc", Frequencies, StatX, StatY,
        StatZ, Impedances, Errors, C);
  }
