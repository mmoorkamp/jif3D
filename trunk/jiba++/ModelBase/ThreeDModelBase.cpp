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

    boost::array<ThreeDModelBase::t3DModelData::index,3> ThreeDModelBase::FindAssociatedIndices(
        const double xcoord, const double ycoord, const double zcoord) const
      {
        const int xindex = std::distance(GetXCoordinates().begin(),
            std::lower_bound(GetXCoordinates().begin(),
                GetXCoordinates().end(), xcoord));
        const int yindex = std::distance(GetYCoordinates().begin(),
            std::lower_bound(GetYCoordinates().begin(),
                GetYCoordinates().end(), ycoord));
        const int zindex = std::distance(GetZCoordinates().begin(),
            std::lower_bound(GetZCoordinates().begin(),
                GetZCoordinates().end(), zcoord));
        boost::array<t3DModelData::index,3> idx =
          {
            { std::max(xindex-1,0), std::max(yindex-1,0), std::max(zindex-1,0) } };
        return idx;
      }

    void ThreeDModelBase::ReadDataFromNetCDF(const NcFile &NetCDFFile,
        const std::string &DataName, const std::string &UnitsName)
      {
        Read3DModelFromNetCDF(NetCDFFile, DataName, UnitsName, XCellSizes,
            YCellSizes, ZCellSizes, Data);
      }

    void ThreeDModelBase::WriteDataToNetCDF(NcFile &NetCDFFile,
        const std::string &DataName, const std::string &UnitsName) const
      {
        Write3DModelToNetCDF(NetCDFFile, DataName, UnitsName, XCellSizes,
            YCellSizes, ZCellSizes, Data);
      }

    void ThreeDModelBase::WriteVTK(std::string filename,
        const std::string &DataName)
      {
        Write3DModelToVTK(filename, DataName, GetXCellSizes(), GetYCellSizes(),
            GetZCellSizes(), GetData());
      }

  }
