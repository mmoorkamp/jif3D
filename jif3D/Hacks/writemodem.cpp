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
#include "../MT/MTData.h"
#include "../MT/TipperData.h"

int main()
  {
    std::string ncfilename = jif3D::AskFilename("Name of MT netcdf file: ");

    jif3D::MTData MTData;
    MTData.ReadNetCDF(ncfilename);
    MTData.WriteModEM(ncfilename+".modem");


    ncfilename = jif3D::AskFilename("Name of Tipper netcdf file: ");

    jif3D::TipperData TipData;
    TipData.ReadNetCDF(ncfilename);
    TipData.WriteModEM(ncfilename+".modem");
  }
