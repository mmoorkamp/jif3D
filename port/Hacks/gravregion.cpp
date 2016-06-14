//============================================================================
// Name        : gravregion.cpp
// Author      : Feb 2, 2011
// Version     : 
// Copyright   : 2011, mmoorkamp
//============================================================================

#include "../Gravity/ReadWriteGravityData.h"
#include "../Global/FileUtil.h"
#include <iostream>

int main()
  {

    jif3D::rvec Data, Error;
    jif3D::ThreeDGravityModel::tMeasPosVec PosX, PosY, PosZ;

    std::string datafilename = jif3D::AskFilename("Data Filename: ");
    //we figure out the type of data (scalar or ftg) from the variables
    //that are in the netcdf file
    jif3D::GravityDataType DataType = jif3D::IdentifyGravityDatafileType(
        datafilename);

    size_t nmeasdata = 1;
    switch (DataType)
      {
    case jif3D::scalar:
      jif3D::ReadScalarGravityMeasurements(datafilename, Data, PosX, PosY, PosZ, Error);
      break;
    case jif3D::ftg:
      nmeasdata = 9;
      jif3D::ReadTensorGravityMeasurements(datafilename, Data, PosX, PosY, PosZ, Error);
      break;
    default:
      //in case we couldn't identify the data in the netcdf file
      //print an error message and exit with an error code
      std::cerr << "Cannot determine the type of data to invert. Aborting."
          << std::endl;
      exit(100);
      break;
      }

    if (Data.empty())
      {
        std::cerr << "No measurements defined" << std::endl;
        exit(100);
      }
    double minx, maxx, miny, maxy;
    std::cout << "Minimum Northing value: ";
    std::cin >> minx;
    std::cout << "Maximum Northing value: ";
    std::cin >> maxx;
    std::cout << "Minimum Easting value: ";
    std::cin >> miny;
    std::cout << "Maximum Easting value: ";
    std::cin >> maxy;

    jif3D::ThreeDGravityModel::tMeasPosVec NewPosX, NewPosY, NewPosZ, NewData;
    const size_t npos = Data.size() / nmeasdata;

    for (size_t i = 0; i < npos; ++i)
      {
        if (PosX.at(i) < maxx && PosX.at(i) > minx && PosY.at(i) < maxy
            && PosY.at(i) > miny)
          {
            NewPosX.push_back(PosX.at(i));
            NewPosY.push_back(PosY.at(i));
            NewPosZ.push_back(PosZ.at(i));
            for (size_t j = 0; j < nmeasdata; ++j)
              {
                NewData.push_back(Data(i * nmeasdata + j));
              }
          }
      }
    jif3D::rvec FinData(NewData.size());
    std::copy(NewData.begin(), NewData.end(), FinData.begin());
    std::string outfilename = jif3D::AskFilename("Output Filename: ", false);
    switch (DataType)
      {
    case jif3D::scalar:
      jif3D::SaveScalarGravityMeasurements(outfilename, FinData, NewPosX,
          NewPosY, NewPosZ, Error);
      break;
    case jif3D::ftg:
      jif3D::SaveTensorGravityMeasurements(outfilename, FinData, NewPosX,
          NewPosY, NewPosZ, Error);
      break;
    default:
      //in case we couldn't identify the data in the netcdf file
      //print an error message and exit with an error code
      std::cerr << "Cannot determine the type of data to invert. Aborting."
          << std::endl;
      exit(100);
      break;
      }
  }
