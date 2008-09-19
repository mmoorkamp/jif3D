//============================================================================
// Name        : ReadWriteGravityData.h
// Author      : Max Moorkamp
// Version     :
// Copyright   : 2008, MM
//============================================================================

#ifndef READWRITEGRAVITYDATA_H_
#define READWRITEGRAVITYDATA_H_

#include "ThreeDGravityModel.h"
#include <string>

namespace jiba
  {
    /** \addtogroup gravity Gravity forward modelling, display and inversion */
     /* @{ */
    //!Read a collection of scalar gravity measurements and associated positions from a netcdf file
    void ReadScalarGravityMeasurements(const std::string &filename,
        ThreeDGravityModel::tScalarMeasVec &Data,
        ThreeDGravityModel::tMeasPosVec &PosX,
        ThreeDGravityModel::tMeasPosVec &PosY,
        ThreeDGravityModel::tMeasPosVec &PosZ);
    //!Save a collection of scalar gravity measurements and associated positions to a netcdf file
    void SaveScalarGravityMeasurements(const std::string &filename,
        const ThreeDGravityModel::tScalarMeasVec &Data,
        const ThreeDGravityModel::tMeasPosVec &PosX,
        const ThreeDGravityModel::tMeasPosVec &PosY,
        const ThreeDGravityModel::tMeasPosVec &PosZ);
    //!Read a collection of tensor gravity measurements and associated positions from a netcdf file
    void ReadTensorGravityMeasurements(const std::string &filename,
        ThreeDGravityModel::tTensorMeasVec &Data,
        ThreeDGravityModel::tMeasPosVec &PosX,
        ThreeDGravityModel::tMeasPosVec &PosY,
        ThreeDGravityModel::tMeasPosVec &PosZ);
    //!Save a collection of tensor gravity measurements and associated positions to a netcdf file
    void SaveTensorGravityMeasurements(const std::string &filename,
        const ThreeDGravityModel::tTensorMeasVec &Data,
        const ThreeDGravityModel::tMeasPosVec &PosX,
        const ThreeDGravityModel::tMeasPosVec &PosY,
        const ThreeDGravityModel::tMeasPosVec &PosZ);
    //! A helper function that reads only measurement positions from a netcdf file regardless of data type
    void ReadMeasPosNetCDF(const std::string filename,
        ThreeDGravityModel::tMeasPosVec &PosX,
        ThreeDGravityModel::tMeasPosVec &PosY,
        ThreeDGravityModel::tMeasPosVec &PosZ);
  /* @} */
  }
#endif /* READWRITEGRAVITYDATA_H_ */
