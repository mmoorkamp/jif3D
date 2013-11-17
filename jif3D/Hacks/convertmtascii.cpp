//============================================================================
// Name        : convertmtt.cpp
// Author      : Sep 21, 2010
// Version     : 
// Copyright   : 2010, mmoorkamp
//============================================================================

#include <iostream>
#include <fstream>
#include <string>
#include "../Global/FileUtil.h"
#include "../Global/FatalException.h"
#include "../MT/ReadWriteImpedances.h"

int main()
{
	std::string filename = jif3D::AskFilename("Ascii file with data: ");
	jif3D::rvec Impedances, Errors;
	std::vector<double> Frequencies, StatXCoord, StatYCoord,
			StatZCoord;

	jif3D::ReadAppResFromAscii(filename, Frequencies, StatXCoord, StatYCoord,
			StatZCoord, Impedances, Errors);
	std::string outfilename = jif3D::AskFilename("Output file: ", false);
	jif3D::WriteImpedancesToNetCDF(outfilename, Frequencies, StatXCoord,
			StatYCoord, StatZCoord, Impedances, Errors);
}
