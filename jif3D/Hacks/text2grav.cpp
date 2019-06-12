//============================================================================
// Name        : gravregion.cpp
// Author      : Feb 2, 2011
// Version     : 
// Copyright   : 2011, mmoorkamp
//============================================================================

#include "../Gravity/ReadWriteGravityData.h"
#include "../Global/FileUtil.h"
#include "../ModelBase/VTKTools.h"
#include <iostream>
#include <fstream>

int main()
  {

    std::vector<double> Data, Error;
    std::vector<double> PosX, PosY, PosZ;

    std::string textfilename = jif3D::AskFilename("Data Filename: ");

    std::ifstream textfile(textfilename.c_str());
    while (textfile.good())
      {
        double currx, curry, currz, currd, curre;
        textfile >> currx >> curry >> currz >> currd >> curre;
        if (textfile.good())
          {
            PosX.push_back(currx);
            PosY.push_back(curry);
            PosZ.push_back(currz);
            Data.push_back(currd);
            Error.push_back(curre);
          }
      }
    std::string outfilename = jif3D::AskFilename("Output Filename: ", false);
    jif3D::SaveScalarGravityMeasurements(outfilename, Data, PosX, PosY, PosZ, Error);
    jif3D::Write3DDataToVTK(outfilename + ".vtk", "grav_accel", Data, PosX, PosY, PosZ);
  }
