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
#include "../MT/MTData.h"

int main()
{
	std::string filename = jif3D::AskFilename("Ascii file with data: ");
	std::vector<double> Impedances, Errors;
	std::vector<double> Frequencies, StatXCoord, StatYCoord,
			StatZCoord;
	jif3D::ReadAppResFromAscii(filename, Frequencies, StatXCoord, StatYCoord,
			StatZCoord, Impedances, Errors);
	jif3D::MTData Data;
	Data.SetMeasurementPoints(StatXCoord,StatYCoord,StatZCoord);
	Data.SetDataAndErrors(Impedances,Errors);
	Data.SetFrequencies(Frequencies);
	Data.CompleteObject();
	std::string outfilename = jif3D::AskFilename("Output file: ", false);
	Data.WriteNetCDF(outfilename);

}
