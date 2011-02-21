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

int main()
    {
    std::string ncfilename = jiba::AskFilename("name of netcdf file; ");

    jiba::rvec Impedances,Errors;
    std::vector<double> Frequencies, StatX, StatY, StatZ;
    jiba::ReadImpedancesFromNetCDF(ncfilename,Frequencies,StatX,StatY,StatZ,Impedances,Errors);
    jiba::WriteImpedancesToMtt(ncfilename,Frequencies,Impedances,Errors);
    }
