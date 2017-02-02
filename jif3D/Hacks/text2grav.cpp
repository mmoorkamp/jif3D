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

    jif3D::rvec Data, Error;
    jif3D::ThreeDGravityModel::tMeasPosVec PosX, PosY, PosZ;
    std::vector<double> DataVec, ErrorVec;

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
            DataVec.push_back(currd);
            ErrorVec.push_back(curre);
          }
      }

    Data.resize(DataVec.size());
    Error.resize(ErrorVec.size());
    std::copy(DataVec.begin(), DataVec.end(), Data.begin());
    std::copy(ErrorVec.begin(), ErrorVec.end(), Error.begin());
    std::string outfilename = jif3D::AskFilename("Output Filename: ", false);
    jif3D::SaveScalarGravityMeasurements(outfilename, Data, PosX, PosY, PosZ, Error);
    jif3D::Write3DDataToVTK(outfilename + ".vtk", "grav_accel", Data, PosX, PosY, PosZ);
  }
