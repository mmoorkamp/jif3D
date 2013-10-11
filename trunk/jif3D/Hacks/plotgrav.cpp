//============================================================================
// Name        : plotgrav.cpp
// Author      : Jul 9, 2010
// Version     : 
// Copyright   : 2010, mmoorkamp
//============================================================================


#include <iostream>
#include <string>
#include "../Gravity/ReadWriteGravityData.h"
#include "../ModelBase/VTKTools.h"
#include "../Global/FileUtil.h"





int main(int argc, char *argv[])
  {

    std::string ScalarFilename = jif3D::AskFilename("Scalar data Filename: ");
    std::string FTGFilename = jif3D::AskFilename("FTG data Filename: ");
    std::vector<double> ScalX, ScalY, ScalZ, FTGX, FTGY, FTGZ;
    jif3D::rvec ScalDat, FTGDat;
    jif3D::ReadScalarGravityMeasurements(ScalarFilename, ScalDat, ScalX, ScalY,
        ScalZ);
    jif3D::ReadTensorGravityMeasurements(FTGFilename, FTGDat, FTGX, FTGY, FTGZ);

    jif3D::Write3DDataToVTK(ScalarFilename + ".vtk", "grav_accel", ScalDat,
        ScalX, ScalY, ScalZ);
    jif3D::Write3DTensorDataToVTK(FTGFilename + ".vtk", "U", FTGDat, FTGX, FTGY,
        FTGZ);

  }
