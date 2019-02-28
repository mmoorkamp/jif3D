//============================================================================
// Name        : ThreeDDCResistivityModel.cpp
// Author      : Zhanjie Shi and Richard.W Hobbs
// Version     : April 2014
// Copyright   : 2014, Zhanjie Shi and Richard.W Hobbs
//============================================================================

#include "ThreeDDCResistivityModel.h"
#include "../Global/FatalException.h"
#include <algorithm>

namespace jif3D
  {
    //we use these names when writing the model to a netcdf file
    static const std::string ResistivityName = "Resistivity";
    static const std::string ResistivityUnit = "ohm.m";

    ThreeDDCResistivityModel::ThreeDDCResistivityModel() :
        SourcePosPosX(), SourcePosPosY(), SourcePosPosZ(), SourceNegPosX(), SourceNegPosY(), SourceNegPosZ(), MeasSecPosX(), MeasSecPosY(), MeasSecPosZ(), SourceIndices()
      {
      }

    ThreeDDCResistivityModel::~ThreeDDCResistivityModel()
      {
      }

    ThreeDDCResistivityModel::ThreeDDCResistivityModel(
        const ThreeDDCResistivityModel &source) :
        ThreeDModelBase(source), SourcePosPosX(source.SourcePosPosX), SourcePosPosY(
            source.SourcePosPosY), SourcePosPosZ(source.SourcePosPosZ), SourceNegPosX(
            source.SourceNegPosX), SourceNegPosY(source.SourceNegPosY), SourceNegPosZ(
            source.SourceNegPosZ), MeasSecPosX(source.MeasSecPosX), MeasSecPosY(
            source.MeasSecPosY), MeasSecPosZ(source.MeasSecPosZ), SourceIndices(
            source.SourceIndices)
      {

      }

    ThreeDDCResistivityModel& ThreeDDCResistivityModel::operator=(const ThreeDModelBase& source)
      {
        if (&source != this)
          {
            ThreeDModelBase::operator=(source);
          }
        return *this;
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
        MeasSecPosX = source.MeasSecPosX;
        MeasSecPosY = source.MeasSecPosY;
        MeasSecPosZ = source.MeasSecPosZ;
        SourceIndices = source.SourceIndices;
        return *this;
      }



    void ThreeDDCResistivityModel::WriteNetCDF(const std::string &filename) const
      {
        netCDF::NcFile DataFile(filename.c_str(), netCDF::NcFile::replace);
        //write the 3D discretized part
        WriteDataToNetCDF(DataFile, ResistivityName, ResistivityUnit);
      }

    void ThreeDDCResistivityModel::ReadNetCDF(const std::string &filename)
      {
        //create the netcdf file object
        netCDF::NcFile DataFile(filename.c_str(), netCDF::NcFile::read);
        //read in the 3D gridded data
        ReadDataFromNetCDF(DataFile, ResistivityName, ResistivityUnit);
      }
  }
