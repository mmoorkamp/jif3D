//============================================================================
// Name        : ModEMModel.cpp
// Author      : 14 Jan 2016
// Version     : 
// Copyright   : 2016, mm489
//============================================================================

#include "ModEMModel.h"

namespace jif3D
  {

    //! Write all model information to a netcdf file
    void ModEMModel::WriteNetCDF(const std::string filename) const
      {
        NcFile DataFile(filename.c_str(), NcFile::Replace);
        //write the 3D discretized part
        WriteDataToNetCDF(DataFile, ConductivityName, ConductivityUnit);
      }
    //! Read all model information from a netcdf file
    void ModEMModel::ReadNetCDF(const std::string filename)
      {
        //create the netcdf file object
        NcFile DataFile(filename.c_str(), NcFile::ReadOnly);
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
      }
    ModEMModel::ModEMModel()
      {
        // TODO Auto-generated constructor stub

      }

    ModEMModel::~ModEMModel()
      {
        // TODO Auto-generated destructor stub
      }

  } /* namespace jif3D */
