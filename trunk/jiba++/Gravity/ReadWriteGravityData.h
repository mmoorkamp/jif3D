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
    void ReadScalarGravityMeasurements(const std::string &filename,
        ThreeDGravityModel::tScalarMeasVec &Data,
        ThreeDGravityModel::tMeasPosVec &PosX,
        ThreeDGravityModel::tMeasPosVec &PosY,
        ThreeDGravityModel::tMeasPosVec &PosZ);

    void SaveScalarGravityMeasurements(const std::string &filename,
        ThreeDGravityModel::tScalarMeasVec &Data,
        ThreeDGravityModel::tMeasPosVec &PosX,
        ThreeDGravityModel::tMeasPosVec &PosY,
        ThreeDGravityModel::tMeasPosVec &PosZ);

    void ReadTensorGravityMeasurements(const std::string &filename,
        ThreeDGravityModel::tTensorMeasVec &Data,
        ThreeDGravityModel::tMeasPosVec &PosX,
        ThreeDGravityModel::tMeasPosVec &PosY,
        ThreeDGravityModel::tMeasPosVec &PosZ);

    void SaveTensorGravityMeasurements(const std::string &filename,
        ThreeDGravityModel::tTensorMeasVec &Data,
        ThreeDGravityModel::tMeasPosVec &PosX,
        ThreeDGravityModel::tMeasPosVec &PosY,
        ThreeDGravityModel::tMeasPosVec &PosZ);

    void ReadMeasPosNetCDF(const std::string filename,
        ThreeDGravityModel::tMeasPosVec &PosX,
        ThreeDGravityModel::tMeasPosVec &PosY,
        ThreeDGravityModel::tMeasPosVec &PosZ);
  }
#endif /* READWRITEGRAVITYDATA_H_ */
