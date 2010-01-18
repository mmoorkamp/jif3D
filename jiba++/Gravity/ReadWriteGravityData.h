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
    /** \addtogroup gravity Gravity forward modeling, display and inversion */
     /* @{ */
    //! This type is used to specify which kind of data a file contains
    enum GravityDataType {scalar, ftg, unknown};
    //! Identify which kind of gravity data is stored in a netcdf file
    GravityDataType IdentifyGravityDatafileType(const std::string &filename);
    //!Read a collection of scalar gravity measurements and associated positions from a netcdf file
    void ReadScalarGravityMeasurements(const std::string &filename,
        jiba::rvec &Data,
        ThreeDGravityModel::tMeasPosVec &PosX,
        ThreeDGravityModel::tMeasPosVec &PosY,
        ThreeDGravityModel::tMeasPosVec &PosZ);
    //!Save a collection of scalar gravity measurements and associated positions to a netcdf file
    void SaveScalarGravityMeasurements(const std::string &filename,
        const jiba::rvec &Data,
        const ThreeDGravityModel::tMeasPosVec &PosX,
        const ThreeDGravityModel::tMeasPosVec &PosY,
        const ThreeDGravityModel::tMeasPosVec &PosZ);
    //!Read a collection of tensor gravity measurements and associated positions from a netcdf file
    void ReadTensorGravityMeasurements(const std::string &filename,
        jiba::rvec &Data,
        ThreeDGravityModel::tMeasPosVec &PosX,
        ThreeDGravityModel::tMeasPosVec &PosY,
        ThreeDGravityModel::tMeasPosVec &PosZ);
    //!Save a collection of tensor gravity measurements and associated positions to a netcdf file
    void SaveTensorGravityMeasurements(const std::string &filename,
        const jiba::rvec &Data,
        const ThreeDGravityModel::tMeasPosVec &PosX,
        const ThreeDGravityModel::tMeasPosVec &PosY,
        const ThreeDGravityModel::tMeasPosVec &PosZ);
  /* @} */
  }
#endif /* READWRITEGRAVITYDATA_H_ */