//============================================================================
// Name        : mtforward.cpp
// Author      : Jul 14, 2009
// Version     : 
// Copyright   : 2009, mmoorkamp
//============================================================================

#include <boost/lambda/lambda.hpp>
#include "../Global/FileUtil.h"
#include "X3DModel.h"
#include "X3DMTCalculator.h"
#include "ReadWriteImpedances.h"

using namespace boost::lambda;
int main()
  {
    std::string modelfilename = jiba::AskFilename("Filename: ", true);
    jiba::X3DModel MTModel;
    MTModel.ReadNetCDF(modelfilename);
    std::string outfilename = jiba::AskFilename("Output filename: ", false);
    jiba::rvec Frequencies(1);
    Frequencies(0) = 10.0;
    MTModel.SetFrequencies().assign(1, 1.0);
    jiba::X3DMTCalculator Calculator;
    jiba::rvec Impedances(Calculator.Calculate(MTModel));
    std::vector<double> XCoord(MTModel.GetConductivities().shape()[0]), YCoord(
        MTModel.GetConductivities().shape()[1]), ZCoord(
        MTModel.GetConductivities().shape()[2]);

    for (size_t i = 0; i < MTModel.GetXCoordinates().size(); ++i)
      {
        for (size_t j = 0; j < MTModel.GetYCoordinates().size(); ++j)
          {
            MTModel.AddMeasurementPoint(MTModel.GetXCoordinates()[i],
                MTModel.GetYCoordinates()[j], 0.0);
          }
      }

    jiba::WriteImpedancesToNetCDF(outfilename, Frequencies,
        MTModel.GetMeasPosX(), MTModel.GetMeasPosY(),
        MTModel.GetMeasPosZ(), Impedances);
  }

