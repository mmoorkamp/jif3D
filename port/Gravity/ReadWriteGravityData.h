//============================================================================
// Name        : ReadWriteGravityData.h
// Author      : Max Moorkamp
// Version     :
// Copyright   : 2008, MM
//============================================================================

#ifndef READWRITEGRAVITYDATA_H_
#define READWRITEGRAVITYDATA_H_

#include <string>

#include "ThreeDGravityModel.h"
#include "../Global/Jif3DGlobal.h"

namespace jif3D
  {
    /** \addtogroup gravity Gravity forward modeling, display and inversion */
     /* @{ */
    //! This type is used to specify which kind of data a file contains
    enum GravityDataType {scalar, ftg, unknown};
    //! Identify which kind of gravity data is stored in a netcdf file
    J3DEXPORT GravityDataType IdentifyGravityDatafileType(const std::string &filename);
    //!Read a collection of scalar gravity measurements and associated positions from a netcdf file
    J3DEXPORT void ReadScalarGravityMeasurements(const std::string &filename,
        jif3D::rvec &Data,
        ThreeDGravityModel::tMeasPosVec &PosX,
        ThreeDGravityModel::tMeasPosVec &PosY,
        ThreeDGravityModel::tMeasPosVec &PosZ,
        jif3D::rvec &Error );
    //!Save a collection of scalar gravity measurements and associated positions to a netcdf file
    J3DEXPORT void SaveScalarGravityMeasurements(const std::string &filename,
        const jif3D::rvec &Data,
        const ThreeDGravityModel::tMeasPosVec &PosX,
        const ThreeDGravityModel::tMeasPosVec &PosY,
        const ThreeDGravityModel::tMeasPosVec &PosZ,
        const jif3D::rvec &Error );
    //!Read a collection of tensor gravity measurements and associated positions from a netcdf file
    J3DEXPORT void ReadTensorGravityMeasurements(const std::string &filename,
        jif3D::rvec &Data,
        ThreeDGravityModel::tMeasPosVec &PosX,
        ThreeDGravityModel::tMeasPosVec &PosY,
        ThreeDGravityModel::tMeasPosVec &PosZ,
        jif3D::rvec &Error);
    //!Save a collection of tensor gravity measurements and associated positions to a netcdf file
    J3DEXPORT void SaveTensorGravityMeasurements(const std::string &filename,
        const jif3D::rvec &Data,
        const ThreeDGravityModel::tMeasPosVec &PosX,
        const ThreeDGravityModel::tMeasPosVec &PosY,
        const ThreeDGravityModel::tMeasPosVec &PosZ,
        const jif3D::rvec &Error );
  /* @} */
  }
#endif /* READWRITEGRAVITYDATA_H_ */
