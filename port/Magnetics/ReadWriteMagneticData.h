//============================================================================
// Name        : ReadWriteMagneticData.h
// Author      : Max Moorkamp
// Version     :
// Copyright   : 2013, MM
//============================================================================

#ifndef READWRITEMAGNETICDATA_H_
#define READWRITEMAGNETICDATA_H_

#include "ThreeDMagneticModel.h"
#include "../Global/VecMat.h"
#include <string>

namespace jif3D
  {
    /** \addtogroup Magnetic Magnetic forward modeling, display and inversion */
     /* @{ */


    //!Read a collection of scalar Magnetic measurements and associated positions from a netcdf file
    void ReadTotalFieldMagneticMeasurements(const std::string &filename,
        jif3D::rvec &Data,
        ThreeDMagneticModel::tMeasPosVec &PosX,
        ThreeDMagneticModel::tMeasPosVec &PosY,
        ThreeDMagneticModel::tMeasPosVec &PosZ,
        jif3D::rvec &Error );
    //!Save a collection of scalar Magnetic measurements and associated positions to a netcdf file
    void SaveTotalFieldMagneticMeasurements(const std::string &filename,
        const jif3D::rvec &Data,
        const ThreeDMagneticModel::tMeasPosVec &PosX,
        const ThreeDMagneticModel::tMeasPosVec &PosY,
        const ThreeDMagneticModel::tMeasPosVec &PosZ,
        const jif3D::rvec &Error );

  /* @} */
  }
#endif /* READWRITEMagneticDATA_H_ */
