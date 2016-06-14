//============================================================================
// Name        : ModEMModel.cpp
// Author      : 14 Jan 2016
// Version     :
// Copyright   : 2016, mm489
//============================================================================

#include "ModEMModel.h"

using netCDF::NcFile;

namespace jif3D
  {

    //! Write all model information to a netcdf file
    void ModEMModel::WriteNetCDF(const std::string &filename) const
      {
        NcFile DataFile(filename, NcFile::replace);
        //write the 3D discretized part
        WriteDataToNetCDF(DataFile, ConductivityName, ConductivityUnit);
      }
    //! Read all model information from a netcdf file
    void ModEMModel::ReadNetCDF(const std::string &filename)
      {
        //create the netcdf file object
        NcFile DataFile(filename, NcFile::read);
        //read in the 3D gridded data
        ReadDataFromNetCDF(DataFile, ConductivityName, ConductivityUnit);
      }
    //! Other models will be copied by the copy operator for the base class
    ModEMModel& ModEMModel::operator=(const ThreeDModelBase& source)
      {
        if (&source != this)
          {
            ThreeDMTModel::operator=(source);
          }

        return *this;
      }

    ModEMModel::ModEMModel()
      {
      }

    ModEMModel::~ModEMModel()
      {
      }

  } /* namespace jif3D */
