//============================================================================
// Name        : ThreeDModelBase.cpp
// Author      : Max Moorkamp
// Version     :
// Copyright   : 2008, MM
//============================================================================

#include "ThreeDModelBase.h"
#include <cassert>
#include <fstream>
#include "NetCDFTools.h"
#include "VTKTools.h"
namespace jiba
  {

    ThreeDModelBase::ThreeDModelBase() :
      XCellSizesChanged(true), YCellSizesChanged(true), ZCellSizesChanged(true)
      {
      }

    ThreeDModelBase::~ThreeDModelBase()
      {
      }


    void ThreeDModelBase::ReadDataFromNetCDF(const NcFile &NetCDFFile,
        const std::string &DataName, const std::string &UnitsName)
      {
        Read3DModelFromNetCDF(NetCDFFile,DataName,UnitsName,XCellSizes,YCellSizes,ZCellSizes,Data);
      }

    void ThreeDModelBase::WriteDataToNetCDF(NcFile &NetCDFFile,
        const std::string &DataName, const std::string &UnitsName) const
      {
        Write3DModelToNetCDF(NetCDFFile,DataName,UnitsName,XCellSizes,YCellSizes,ZCellSizes,Data);
      }



    void ThreeDModelBase::WriteVTK(std::string filename,
        const std::string &DataName)
      {
        Write3DModelToVTK(filename,DataName,GetXCellSizes(),GetYCellSizes(),GetZCellSizes(),GetData());
      }

  }
