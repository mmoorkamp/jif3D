//============================================================================
// Name        : ReadWriteSparseMatrix.h
// Author      : 1 Jul 2013
// Version     : 
// Copyright   : 2013, mm489
//============================================================================

#include <string>
#include "VecMat.h"

#ifndef READWRITECOVARIANCE_H_
#define READWRITECOVARIANCE_H_

namespace jif3D
  {
    void ReadSparseMatrixFromNetcdf(const std::string &filename, jif3D::comp_mat &CoVar, const std::string &Name);
    void WriteCovarianceToNetcdf(const std::string &filename,
        const jif3D::comp_mat &CoVar, const std::string &Name);

  }

#endif /*READWRITECOVARIANCE_H_*/
