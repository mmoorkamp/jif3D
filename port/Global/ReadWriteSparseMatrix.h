//============================================================================
// Name        : ReadWriteSparseMatrix.h
// Author      : 1 Jul 2013
// Version     : 
// Copyright   : 2013, mm489
//============================================================================

#ifndef READWRITECOVARIANCE_H_
#define READWRITECOVARIANCE_H_

#include <string>

#include "VecMat.h"
#include "Jif3DGlobal.h"

namespace jif3D
  {
    J3DEXPORT void ReadSparseMatrixFromNetcdf(const std::string &filename, jif3D::comp_mat &CoVar, const std::string &Name);
    J3DEXPORT void WriteSparseMatrixToNetcdf(const std::string &filename,
        const jif3D::comp_mat &CoVar, const std::string &Name);
  }

#endif /*READWRITECOVARIANCE_H_*/
