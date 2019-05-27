//============================================================================
// Name        : ReadWriteGravityData.h
// Author      : Max Moorkamp
// Version     :
// Copyright   : 2008, MM
//============================================================================

#ifndef READWRITEGRAVITYDATA_H_
#define READWRITEGRAVITYDATA_H_

#include "../Global/Jif3DGlobal.h"
#include "../Global/VecMat.h"
#include <string>
#include <vector>

namespace jif3D
  {
    /** \addtogroup gravity Gravity forward modeling, display and inversion */
    /* @{ */
    //! This type is used to specify which kind of data a file contains
    enum GravityDataType
      {
      scalar, ftg, unknown
      };
    //! Identify which kind of gravity data is stored in a netcdf file
    J3DEXPORT GravityDataType IdentifyGravityDatafileType(const std::string &filename);
    //!Read a collection of scalar gravity measurements and associated positions from a netcdf file
    J3DEXPORT void ReadScalarGravityMeasurements(const std::string &filename,
        std::vector<double>  &Data, std::vector<double> &PosX, std::vector<double> &PosY,
        std::vector<double> &PosZ, std::vector<double>  &Error);
    //!Save a collection of scalar gravity measurements and associated positions to a netcdf file
    J3DEXPORT void SaveScalarGravityMeasurements(const std::string &filename,
        const std::vector<double>  &Data, const std::vector<double> &PosX,
        const std::vector<double> &PosY, const std::vector<double> &PosZ,
        const std::vector<double>  &Error);
    //!Read a collection of tensor gravity measurements and associated positions from a netcdf file
    J3DEXPORT void ReadTensorGravityMeasurements(const std::string &filename,
        std::vector<double>  &Data, std::vector<double> &PosX, std::vector<double> &PosY,
        std::vector<double> &PosZ, std::vector<double>  &Error);
    //!Save a collection of tensor gravity measurements and associated positions to a netcdf file
    J3DEXPORT void SaveTensorGravityMeasurements(const std::string &filename,
        const std::vector<double>  &Data, const std::vector<double> &PosX,
        const std::vector<double> &PosY, const std::vector<double> &PosZ,
        const std::vector<double>  &Error);
  /* @} */
  }
#endif /* READWRITEGRAVITYDATA_H_ */
