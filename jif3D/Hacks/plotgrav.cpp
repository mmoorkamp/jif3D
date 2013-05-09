//============================================================================
// Name        : plotgrav.cpp
// Author      : Jul 9, 2010
// Version     : 
// Copyright   : 2010, mmoorkamp
//============================================================================


#include <iostream>
#include <string>
#include "ReadWriteGravityData.h"
#include "../ModelBase/VTKTools.h"
#include "../Global/FileUtil.h"





int main(int argc, char *argv[])
  {

    std::string ScalarFilename = jiba::AskFilename("Scalar data Filename: ");
    std::string FTGFilename = jiba::AskFilename("FTG data Filename: ");
    std::vector<double> ScalX, ScalY, ScalZ, FTGX, FTGY, FTGZ;
    jiba::rvec ScalDat, FTGDat;
    jiba::ReadScalarGravityMeasurements(ScalarFilename, ScalDat, ScalX, ScalY,
        ScalZ);
    jiba::ReadTensorGravityMeasurements(FTGFilename, FTGDat, FTGX, FTGY, FTGZ);

    jiba::Write3DDataToVTK(ScalarFilename + ".vtk", "grav_accel", ScalDat,
        ScalX, ScalY, ScalZ);
    jiba::Write3DTensorDataToVTK(FTGFilename + ".vtk", "U", FTGDat, FTGX, FTGY,
        FTGZ);

  }
