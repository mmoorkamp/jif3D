//============================================================================
// Name        : ThreeDDCResistivityModel.cpp
// Author      : Zhanjie Shi and Richard.W Hobbs
// Version     : April 2014
// Copyright   : 2014, Zhanjie Shi and Richard.W Hobbs
//============================================================================

#include "ThreeDDCResistivityModel.h"
#include "../Global/FatalException.h"
#include <boost/bind.hpp>
#include <boost/numeric/conversion/cast.hpp>
#include <cassert>
#include <algorithm>

namespace jif3D
  {
    //we use these names when writing the model to a netcdf file
    static const std::string ResistivityName = "Resistivity";
    static const std::string ResistivityUnit = "ohm.m";

    ThreeDDCResistivityModel::ThreeDDCResistivityModel() :
      SourcePosPosX(), SourcePosPosY(), SourcePosPosZ(), SourceNegPosX(), SourceNegPosY(), SourceNegPosZ(),
          SourceIndices(), MeasSecPosX(), MeasSecPosY(), MeasSecPosZ()
      {
      }

    ThreeDDCResistivityModel::~ThreeDDCResistivityModel()
      {
      }

    ThreeDDCResistivityModel::ThreeDDCResistivityModel(const ThreeDDCResistivityModel &source) :
      ThreeDModelBase(source), SourcePosPosX(source.SourcePosPosX), SourcePosPosY(source.SourcePosPosY),
      SourcePosPosZ(source.SourcePosPosZ), SourceNegPosX(source.SourceNegPosX), SourceNegPosY(source.SourceNegPosY),
      SourceNegPosZ(source.SourceNegPosZ), SourceIndices(source.SourceIndices)
      {

      }

    ThreeDDCResistivityModel& ThreeDDCResistivityModel::operator=(
        const ThreeDDCResistivityModel& source)
      {
        if (this == &source)
          return *this;
        ThreeDModelBase::operator =(source);
        SourcePosPosX = source.SourcePosPosX;
        SourcePosPosY = source.SourcePosPosY;
        SourcePosPosZ = source.SourcePosPosZ;
        SourceNegPosX = source.SourceNegPosX;
        SourceNegPosY = source.SourceNegPosY;
        SourceNegPosZ = source.SourceNegPosZ;
        SourceIndices = source.SourceIndices;
        return *this;
      }

    void ThreeDDCResistivityModel::SetOrigin(const double x, const double y,
        const double z)
      {
        //transform the source coordinates from old model to real coordinates
        //the coordinates of the receivers are changed by the implementation
        //in the base class that we call below
        std::transform(SourcePosPosX.begin(), SourcePosPosX.end(),
            SourcePosPosX.begin(), boost::bind(std::plus<double>(), _1, XOrigin
                - x));
        std::transform(SourcePosPosY.begin(), SourcePosPosY.end(),
            SourcePosPosY.begin(), boost::bind(std::plus<double>(), _1, YOrigin
                - y));
        std::transform(SourcePosPosZ.begin(), SourcePosPosZ.end(),
            SourcePosPosZ.begin(), boost::bind(std::plus<double>(), _1, ZOrigin
                - z));
        std::transform(SourceNegPosX.begin(), SourceNegPosX.end(),
            SourceNegPosX.begin(), boost::bind(std::plus<double>(), _1, XOrigin
                - x));
        std::transform(SourceNegPosY.begin(), SourceNegPosY.end(),
            SourceNegPosY.begin(), boost::bind(std::plus<double>(), _1, YOrigin
                - y));
        std::transform(SourceNegPosZ.begin(), SourceNegPosZ.end(),
            SourceNegPosZ.begin(), boost::bind(std::plus<double>(), _1, ZOrigin
                - z));
        std::transform(MeasSecPosX.begin(), MeasSecPosX.end(),
        		MeasSecPosX.begin(), boost::bind(std::plus<double>(), _1, XOrigin
                - x));
        std::transform(MeasSecPosY.begin(), MeasSecPosY.end(),
        		MeasSecPosY.begin(), boost::bind(std::plus<double>(), _1, YOrigin
                - y));
        std::transform(MeasSecPosZ.begin(), MeasSecPosZ.end(),
        		MeasSecPosZ.begin(), boost::bind(std::plus<double>(), _1, ZOrigin
                - z));
        //we have to call the base implementation in the end because
        //it changes the first measurement positions and the Origin
        ThreeDModelBase::SetOrigin(x, y, z);
      }

    void ThreeDDCResistivityModel::WriteNetCDF(const std::string filename) const
      {
        NcFile DataFile(filename.c_str(), NcFile::Replace);
        //write the 3D discretized part
        WriteDataToNetCDF(DataFile, ResistivityName, ResistivityUnit);
      }

    void ThreeDDCResistivityModel::ReadNetCDF(const std::string filename)
      {
        //create the netcdf file object
        NcFile DataFile(filename.c_str(), NcFile::ReadOnly);
        //read in the 3D gridded data
        ReadDataFromNetCDF(DataFile, ResistivityName, ResistivityUnit);
      }
  }
