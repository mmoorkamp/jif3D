//============================================================================
// Name        : NetCDFPortHelper.cpp
// Author      : Ashley Davis (ad299) https://github.com/SgtCoDFish
// Version     :
// Copyright   : 2016, AD
//============================================================================

#include <vector>

#include "NetCDFPortHelper.h"

using netCDF::NcVar;
using netCDF::NcDim;

namespace jif3D
{
namespace cxxport
{

// Implemented in header file now to keep helper header-only
//long num_vals(const NcVar &var)
//{
//  const std::vector<NcDim> dims = var.getDims();
//
//  long prod = 1;
//
//  for (size_t i = 0u; i < dims.size(); ++i)
//  {
//    prod *= dims[i].getSize();
//  }
//
//  return prod;
//}

}
}
