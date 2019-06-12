//============================================================================
// Name        : gravregion.cpp
// Author      : Feb 2, 2011
// Version     : 
// Copyright   : 2011, mmoorkamp
//============================================================================

#include "../Gravity/ReadWriteGravityData.h"
#include "../Gravity/ScalarGravityData.h"
#include "../Gravity/TensorGravityData.h"
#include "../Global/FileUtil.h"
#include <boost/shared_ptr.hpp>
#include <iostream>

int main()
  {

    boost::shared_ptr<jif3D::GeneralData> Data;

    std::string datafilename = jif3D::AskFilename("Data Filename: ");
    //we figure out the type of data (scalar or ftg) from the variables
    //that are in the netcdf file
    jif3D::GravityDataType DataType = jif3D::IdentifyGravityDatafileType(datafilename);

    size_t nmeasdata = 1;
    switch (DataType)
      {
    case jif3D::scalar:
      Data = boost::make_shared<jif3D::ScalarGravityData>();
      break;
    case jif3D::ftg:
      nmeasdata = 9;
      Data = boost::make_shared<jif3D::TensorGravityData>();
      break;
    default:
      //in case we couldn't identify the data in the netcdf file
      //print an error message and exit with an error code
      std::cerr << "Cannot determine the type of data to invert. Aborting." << std::endl;
      exit(100);
      break;
      }
    Data->ReadNetCDF(datafilename);

    if (Data->GetData().empty())
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

    std::vector<double> NewPosX, NewPosY, NewPosZ, NewData, NewError;
    const size_t npos = Data->GetData().size() / nmeasdata;

    for (size_t i = 0; i < npos; ++i)
      {
        if (Data->GetMeasPosX().at(i) < maxx && Data->GetMeasPosX().at(i) > minx && Data->GetMeasPosY().at(i) < maxy
            && Data->GetMeasPosY().at(i) > miny)
          {
            NewPosX.push_back(Data->GetMeasPosX().at(i));
            NewPosY.push_back(Data->GetMeasPosY().at(i));
            NewPosZ.push_back(Data->GetMeasPosZ().at(i));
            for (size_t j = 0; j < nmeasdata; ++j)
              {
                NewData.push_back(Data->GetData().at(i * nmeasdata + j));
                NewError.push_back(Data->GetErrors().at(i * nmeasdata + j));
              }
          }
      }
    Data->ClearMeasurementPoints();
    Data->SetMeasurementPoints(NewPosX,NewPosY,NewPosZ);
    Data->SetDataAndErrors(NewData,NewError);
    Data->WriteNetCDF(datafilename+".rot.nc");
    Data->WriteMeasurementPoints(datafilename+".rot.vtk");
  }
