//============================================================================
// Name        : ThreeDModelBase.cpp
// Author      : Max Moorkamp
// Version     :
// Copyright   : 2008, MM
//============================================================================

#include "ThreeDModelBase.h"
#include <cassert>
#include <fstream>
#include <boost/bind.hpp>
#include "NetCDFTools.h"
#include "VTKTools.h"
namespace jiba
  {

    ThreeDModelBase::ThreeDModelBase() :
      XOrigin(0.0), YOrigin(0.0), ZOrigin(0.0), XCellSizesChanged(true),
          YCellSizesChanged(true), ZCellSizesChanged(true)
      {
      }

    ThreeDModelBase::~ThreeDModelBase()
      {
      }

    boost::array<ThreeDModelBase::t3DModelData::index, 3> ThreeDModelBase::FindAssociatedIndices(
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
        boost::array<t3DModelData::index, 3> idx =
          {
            { std::max(xindex - 1, 0), std::max(yindex - 1, 0), std::max(zindex
                - 1, 0) } };
        return idx;
      }

    void ThreeDModelBase::SetOrigin(const double x, const double y,
        const double z)
      {
        //transform the measurement coordinates from old model to real coordinates
        std::transform(MeasPosX.begin(), MeasPosX.end(), MeasPosX.begin(),
            boost::bind(std::plus<double>(), _1, XOrigin));
        std::transform(MeasPosY.begin(), MeasPosY.end(), MeasPosY.begin(),
            boost::bind(std::plus<double>(), _1, YOrigin));
        std::transform(MeasPosZ.begin(), MeasPosZ.end(), MeasPosZ.begin(),
            boost::bind(std::plus<double>(), _1, ZOrigin));
        //now transform from real to new model coordinates
        std::transform(MeasPosX.begin(), MeasPosX.end(), MeasPosX.begin(),
            boost::bind(std::minus<double>(), _1, x));
        std::transform(MeasPosY.begin(), MeasPosY.end(), MeasPosY.begin(),
            boost::bind(std::minus<double>(), _1, y));
        std::transform(MeasPosZ.begin(), MeasPosZ.end(), MeasPosZ.begin(),
            boost::bind(std::minus<double>(), _1, z));
        XOrigin = x;
        YOrigin = y;
        ZOrigin = z;

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

    void ThreeDModelBase::ReadMeasPosNetCDF(const std::string filename)
      {
        jiba::ReadMeasPosNetCDF(filename, MeasPosX, MeasPosY, MeasPosZ);
      }

    void ThreeDModelBase::ReadMeasPosAscii(const std::string filename)
      {
        std::ifstream infile(filename.c_str());
        double posx, posy, posz;
        while (infile.good())
          {
            infile >> posx >> posy >> posz;
            if (infile.good())
              {
                MeasPosX.push_back(posx);
                MeasPosY.push_back(posy);
                MeasPosZ.push_back(posz);
              }
          }
        assert(MeasPosX.size() == MeasPosY.size());
        assert(MeasPosX.size() == MeasPosZ.size());
      }

  }
