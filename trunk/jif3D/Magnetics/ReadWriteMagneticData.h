//============================================================================
// Name        : ReadWriteMagneticData.h
// Author      : Max Moorkamp
// Version     :
// Copyright   : 2013, MM
//============================================================================

#ifndef READWRITEMAGNETICDATA_H_
#define READWRITEMAGNETICDATA_H_

#include <string>
#include <vector>

#include "../Global/VecMat.h"
#include "../Global/Jif3DGlobal.h"
#include "ThreeDSusceptibilityModel.h"

namespace jif3D
  {
    /** \addtogroup Magnetic Magnetic forward modeling, display and inversion */
    /* @{ */

    //!Read a collection of scalar Magnetic measurements and associated positions from a netcdf file
    J3DEXPORT void ReadTotalFieldMagneticMeasurements(const std::string &filename,
        std::vector<double> &Data, std::vector<double> &PosX, std::vector<double> &PosY,
        std::vector<double> &PosZ, std::vector<double> &Error);
    //!Save a collection of scalar Magnetic measurements and associated positions to a netcdf file
    J3DEXPORT void SaveTotalFieldMagneticMeasurements(const std::string &filename,
        const std::vector<double> &Data, const std::vector<double> &PosX,
        const std::vector<double> &PosY, const std::vector<double> &PosZ,
        const std::vector<double> &Error);
    //! Read three component magnetic data
    void ReadMagneticComponentMeasurements(const std::string &filename, jif3D::rvec &Data,
        std::vector<double> &PosX, std::vector<double> &PosY, std::vector<double> &PosZ,
        jif3D::rvec &Error);
    //! Write three component magnetic data
    void SaveMagneticComponentMeasurements(const std::string &filename,
        const jif3D::rvec &Data, const std::vector<double> &PosX,
        const std::vector<double> &PosY, const std::vector<double> &PosZ,
        const jif3D::rvec &Error);

  /* @} */
  }
#endif /* READWRITEMagneticDATA_H_ */
