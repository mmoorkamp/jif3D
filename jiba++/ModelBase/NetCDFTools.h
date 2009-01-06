//============================================================================
// Name        : NetCDFTools.h
// Author      : Max Moorkamp
// Version     :
// Copyright   : 2008, MM
//============================================================================

#ifndef NETCDFTOOLS_H_
#define NETCDFTOOLS_H_

#include <string>
#include <netcdfcpp.h>
#include "ThreeDModelBase.h"

namespace jiba
  {
    /** \addtogroup modelbase Basic classes and routines for 3D models */
    /* @{ */

    //! Read a single CellSize variable from a netcdf files
    size_t ReadSizesFromNetCDF(const NcFile &NetCDFFile,
        const std::string &SizeName, ThreeDModelBase::t3DModelDim &CellSize);
    //! Write the length of the model cell along a single dimension to the file
    NcDim *WriteSizesToNetCDF(NcFile &NetCDFFile, const std::string &SizeName,
        const ThreeDModelBase::t3DModelDim &CellSize);
    //! Read a 3D model from netcdf, this is the preferred storage format
    void Read3DModelFromNetCDF(const NcFile &NetCDFFile,
        const std::string &DataName, const std::string &UnitsName,
        ThreeDModelBase::t3DModelDim &XCellSizes,
        ThreeDModelBase::t3DModelDim &YCellSizes,
        ThreeDModelBase::t3DModelDim &ZCellSizes,
        ThreeDModelBase::t3DModelData &Data);
    //! Write a 3D model to a netcdf file, this is the preferred storage format
    void Write3DModelToNetCDF(NcFile &NetCDFFile, const std::string &DataName,
        const std::string &UnitsName,
        const ThreeDModelBase::t3DModelDim &XCellSizes,
        const ThreeDModelBase::t3DModelDim &YCellSizes,
        const ThreeDModelBase::t3DModelDim &ZCellSizes,
        const ThreeDModelBase::t3DModelData &Data);
    /* @} */
  }
#endif /*NETCDFTOOLS_H_*/
