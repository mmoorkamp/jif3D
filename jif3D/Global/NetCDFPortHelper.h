//============================================================================
// Name        : NetCDFPortHelper.h
// Author      : Ashley Davis (ad299) https://github.com/SgtCoDFish
// Version     :
// Copyright   : 2016, AD
//============================================================================

/*
 * This file contains various functions aimed to simplify the port from the legacy
 * netcdf-cxx package to the newer version.
 *
 * Some methods that were provided by the old API are deprecated or otherwise
 * not supported; this file largely recreates them to minimise the need to
 * rewrite huge swathes of code.
 *
 * NOTE THAT NEWLY WRITTEN CODE SHOULD LIKELY NOT USE THIS FILE:
 * There are new, better ways to write most of these functions. This file is purely
 * to make porting more smooth.
 */

#ifndef GLOBAL_NETCDFPORTHELPER_H_
#define GLOBAL_NETCDFPORTHELPER_H_

#include <netcdf>
#include <vector>

#include "Jif3DGlobal.h"

namespace jif3D
{
namespace cxxport
{

/*! Analogous to the old "long NcVar::num_vals(void) const" method.
 *
 * This function returns the product of the sizes of the dimensions of the given variable.
 * @param var the initialised NcVar from which the dimension sizes will be obtained.
 */
template<typename T = long> J3DEXPORT T num_vals(const netCDF::NcVar &var)
{
  const std::vector<netCDF::NcDim> dims = var.getDims();

  T prod = 1;

  for (size_t i = 0u; i < dims.size(); ++i)
  {
    prod *= dims[i].getSize();
  }

  return prod;
}

/*!
 * Analogous to old bool NcVar::put(T * vals, long c0, long c1, ...) method.
 *
 * Original documentation as follows:
 *
 * "Put scalar or 1, ..., 5 dimensional arrays by providing enough
 * arguments.  Arguments are edge lengths, and their number must not
 * exceed variable's dimensionality.  Start corner is [0,0,..., 0] by
 * default, but may be reset using the set_cur() member.  FALSE is
 * returned if type of values does not match type for variable."
 *
 * set_cur() is not relevant here. Its functionality is described below:
 * Upon creation of a legacy NcVar object, the_cur is initialised to a long[]
 * with a size of NC_MAX_DIMS (==1024 in at least some recent versions of NetCDF,
 * defined in the C library's netcdf.h)
 *
 * This array is initialised to all 0s, and can be changed or reset as needed by
 * NcVar::set_cur(long * cur) or NcVar::set_cur(long, long, long, long, long).
 * These methods are only obviously used in the "put_rec" method which doesn't
 * seem to be called in Jif3D.
 *
 * As such, the cursor can basically be reduced to an long-type array of 0s, and
 * is treated as such in the NetCDFPortHelper.
 */
template<typename T> J3DEXPORT bool put_legacy_ncvar(netCDF::NcVar &var,
    const T * const buffer, const long c1 = 0L, const long c2 = 0L,
    const long c3 = 0L, const long c4 = 0L, const long c5 = 0L)
{
  std::vector<size_t> start(5, 0UL);
  std::vector<size_t> count;

  if (c1 > 0L)
  {
    count.push_back(c1);
  }

  if (c2 > 0L)
  {
    count.push_back(c2);
  }

  if (c3 > 0L)
  {
    count.push_back(c3);
  }

  if (c4 > 0L)
  {
    count.push_back(c4);
  }

  if (c5 > 0L)
  {
    count.push_back(c5);
  }

  var.putVar(start, count, buffer);
  return true;
}

/*!
 * Analogous to old bool NcVar::get(T * vals, long c0, long c1, ...) const method.
 *
 * Original documentation was as follows:
 *
 * "Get scalar or 1, ..., 5 dimensional arrays by providing enough
 * arguments.  Arguments are edge lengths, and their number must not
 * exceed variable's dimensionality.  Start corner is [0,0,..., 0] by
 * default, but may be reset using the set_cur() member."
 *
 * (set_cur() is not relevant here, see note on put_legacy_ncvar)
 */
template<typename T> J3DEXPORT bool get_legacy_ncvar(const netCDF::NcVar &var, T * buffer,
    const unsigned long c1 = 0L, const unsigned long c2 = 0L,
    const unsigned long c3 = 0L, const unsigned long c4 = 0L,
    const unsigned long c5 = 0L)
{
  std::vector<size_t> start(5, 0UL);
  std::vector<size_t> count;

  if (c1 > 0L)
  {
    count.push_back(c1);
  }

  if (c2 > 0L)
  {
    count.push_back(c2);
  }

  if (c3 > 0L)
  {
    count.push_back(c3);
  }

  if (c4 > 0L)
  {
    count.push_back(c4);
  }

  if (c5 > 0L)
  {
    count.push_back(c5);
  }

  var.getVar(start, count, buffer);
  return true;
}

/*!
 * Analogous to the old "long * NcVar::edges() const" method.
 * This version changes the return type to allow RAII and not rely
 * on the programmer to remember to delete[] the return like the original
 * did.
 */
template<typename T = long> J3DEXPORT std::vector<T> get_legacy_var_edges(
    const netCDF::NcVar &var)
{
  const std::vector<netCDF::NcDim> dims = var.getDims();
  std::vector<T> ret;
  ret.reserve(dims.size());

  for(size_t i = 0u; i < dims.size(); i++) {
      ret.push_back(dims[i].getSize());
  }

  return ret;
}

}

}

#endif
