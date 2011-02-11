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

    jiba::rvec Data;
    jiba::ThreeDGravityModel::tMeasPosVec PosX, PosY, PosZ;

    std::string datafilename = jiba::AskFilename("Data Filename: ");
    //we figure out the type of data (scalar or ftg) from the variables
    //that are in the netcdf file
    jiba::GravityDataType DataType = jiba::IdentifyGravityDatafileType(
        datafilename);

    size_t nmeasdata = 1;
    switch (DataType)
      {
    case jiba::scalar:
      jiba::ReadScalarGravityMeasurements(datafilename, Data, PosX, PosY, PosZ);
      break;
    case jiba::ftg:
      nmeasdata = 9;
      jiba::ReadTensorGravityMeasurements(datafilename, Data, PosX, PosY, PosZ);
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

    jiba::ThreeDGravityModel::tMeasPosVec NewPosX, NewPosY, NewPosZ, NewData;
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
    jiba::rvec FinData(NewData.size());
    std::copy(NewData.begin(), NewData.end(), FinData.begin());
    std::string outfilename = jiba::AskFilename("Output Filename: ", false);
    switch (DataType)
      {
    case jiba::scalar:
      jiba::SaveScalarGravityMeasurements(outfilename, FinData, NewPosX,
          NewPosY, NewPosZ);
      break;
    case jiba::ftg:
      jiba::SaveTensorGravityMeasurements(outfilename, FinData, NewPosX,
          NewPosY, NewPosZ);
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
