//============================================================================
// Name        : plotgrav.cpp
// Author      : Jul 9, 2010
// Version     : 
// Copyright   : 2010, mmoorkamp
//============================================================================


#include <iostream>
#include <string>
#include "../Gravity/ScalarGravityData.h"
#include "../Gravity/TensorGravityData.h"
#include "../ModelBase/VTKTools.h"
#include "../Global/FileUtil.h"





int main(int argc, char *argv[])
  {

    std::string ScalarFilename = jif3D::AskFilename("Scalar data Filename: ");
    std::string FTGFilename = jif3D::AskFilename("FTG data Filename: ");
    jif3D::ScalarGravityData ScalData;
    jif3D::TensorGravityData TensData;

    ScalData.ReadNetCDF(ScalarFilename);
    ScalData.WriteVTK(ScalarFilename + ".vtk");
    TensData.ReadNetCDF(FTGFilename);
    TensData.WriteVTK(FTGFilename + ".vtk");

  }
