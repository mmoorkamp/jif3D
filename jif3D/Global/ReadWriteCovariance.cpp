//============================================================================
// Name        : ReadWriteCovariance.cpp
// Author      : 1 Jul 2013
// Version     : 
// Copyright   : 2013, mm489
//============================================================================

#include <netcdfcpp.h>
#include "ReadWriteCovariance.h"

namespace jif3D
  {
    const std::string CovarName = "Covariance";

    void ReadCovarianceFromNetcdf(const std::string &filename,
        const jif3D::comp_mat &CoVar)
      {
        NcFile DataFile(filename.c_str(), NcFile::ReadOnly);
      }

    void WriteCovarianceToNetcdf(const std::string &filename, jif3D::comp_mat &CoVar)
      {
        NcFile DataFile(filename.c_str(), NcFile::Replace);
      }

  }
