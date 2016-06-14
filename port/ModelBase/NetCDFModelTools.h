//============================================================================
// Name        : NetCDFModelTools.h
// Author      : Max Moorkamp
// Version     :
// Copyright   : 2008, MM
//============================================================================

#ifndef NetCDFModelTools_H_
#define NetCDFModelTools_H_

#include <string>
#include <netcdf>
#include "../Global/Jif3DGlobal.h"
#include "ThreeDModelBase.h"

/*! \file NetCDFModelTools.h
 * This file contains various functions to read and write 3D gridded models to and from
 * netcdf files. These functions are used by ThreeDModelBase and derived classes, but
 * are also useful on their own.
 */

namespace jif3D
  {
    /** \addtogroup modelbase Basic classes and routines for 3D models */
    /* @{ */
    //! The name used for the x-coordinates of the measurements in netcdf files
    static const std::string MeasPosXName = "MeasPosX";
    //! The name used for the y-coordinates of the measurements in netcdf files
    static const std::string MeasPosYName = "MeasPosY";
    //! The name used for the z-coordinates of the measurements in netcdf files
    static const std::string MeasPosZName = "MeasPosZ";
    //! The name used for the index of the measurements in netcdf files
    static const std::string StationNumberName = "StationNumber";
    //! Read a 3D model from netcdf, this is the preferred storage format
    J3DEXPORT void Read3DModelFromNetCDF(const netCDF::NcFile &NetCDFFile,
        const std::string &DataName, const std::string &UnitsName,
        ThreeDModelBase::t3DModelDim &XCellSizes,
        ThreeDModelBase::t3DModelDim &YCellSizes,
        ThreeDModelBase::t3DModelDim &ZCellSizes,
        ThreeDModelBase::t3DModelData &Data);
    //! Write a 3D model to a netcdf file using an opened NcFile, this is the preferred storage format
    J3DEXPORT void Write3DModelToNetCDF(netCDF::NcFile &NetCDFFile, const std::string &DataName,
        const std::string &UnitsName,
        const ThreeDModelBase::t3DModelDim &XCellSizes,
        const ThreeDModelBase::t3DModelDim &YCellSizes,
        const ThreeDModelBase::t3DModelDim &ZCellSizes,
        const ThreeDModelBase::t3DModelData &Data);
    //! Write a 3D model to a netcdf file, this version creates the necessary NcFile object
    inline void Write3DModelToNetCDF(const std::string filename,
        const std::string &DataName, const std::string &UnitsName,
        const ThreeDModelBase::t3DModelDim &XCellSizes,
        const ThreeDModelBase::t3DModelDim &YCellSizes,
        const ThreeDModelBase::t3DModelDim &ZCellSizes,
        const ThreeDModelBase::t3DModelData &Data)
      {
      netCDF::NcFile DataFile(filename, netCDF::NcFile::replace);
        Write3DModelToNetCDF(DataFile, DataName, UnitsName, XCellSizes,
            YCellSizes, ZCellSizes, Data);
      }
    //! A helper function that reads only measurement positions from a netcdf file regardless of data type
    J3DEXPORT void ReadMeasPosNetCDF(const std::string &filename,
        ThreeDModelBase::tMeasPosVec &PosX, ThreeDModelBase::tMeasPosVec &PosY,
        ThreeDModelBase::tMeasPosVec &PosZ);

  /* @} */
  }
#endif /*NetCDFModelTools_H_*/
