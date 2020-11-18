/*
 * ThreeDDCResistivityModel.cpp
 *
 *  Created on: Nov 13, 2020
 *      Author: zhanjie
 */

#include "../DCResistivity/ThreeDDCResistivityModel.h"

#include "../Global/FatalException.h"
#include <boost/numeric/conversion/cast.hpp>
#include <cassert>
#include <algorithm>
#include <netcdf>

using netCDF::NcFile; // @suppress("Symbol is not resolved")
using netCDF::NcVar; // @suppress("Symbol is not resolved")
using netCDF::NcDim; // @suppress("Symbol is not resolved")

namespace jif3D
  {
    //we use these names when writing the model to a netcdf file
    static const std::string ResistivityName = "Resistivity";
    static const std::string ResistivityUnit = "ohm.m";

    ThreeDDCResistivityModel::ThreeDDCResistivityModel()
      {
      }

    ThreeDDCResistivityModel::~ThreeDDCResistivityModel()
      {
      }

    ThreeDDCResistivityModel& ThreeDDCResistivityModel::operator=(const ThreeDModelBase& source)
      {
        if (this == &source)
          return *this;
        ThreeDModelBase::operator =(source);

        return *this;
      }


    void ThreeDDCResistivityModel::WriteNetCDF(const std::string &filename) const
      {
        NcFile DataFile(filename.c_str(), NcFile::replace); // @suppress("Type cannot be resolved") // @suppress("Symbol is not resolved")
        //write the 3D discretized part
        WriteDataToNetCDF(DataFile, ResistivityName, ResistivityUnit); // @suppress("Invalid arguments")
      }

    void ThreeDDCResistivityModel::ReadNetCDF(const std::string &filename)
      {
        //create the netcdf file object
        NcFile DataFile(filename.c_str(), NcFile::read); // @suppress("Type cannot be resolved") // @suppress("Symbol is not resolved")
        //read in the 3D gridded data
        ReadDataFromNetCDF(DataFile, ResistivityName, ResistivityUnit); // @suppress("Invalid arguments")
      }
  }
