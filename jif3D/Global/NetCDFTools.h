//============================================================================
// Name        : NetCDFTools.h
// Author      : 1 Jul 2013
// Version     :
// Copyright   : 2013, mm489
//============================================================================

#ifndef NETCDFTOOLS_H_
#define NETCDFTOOLS_H_

#include <netcdf>
#include <string>

#include "NetCDFPortHelper.h"
#include "Jif3DGlobal.h"

namespace jif3D
  {

    //! Read a vector from a netcdf file
    /*! This helper function reads a single vectorial quantity from
     *  an already opened netcdf file. It is used to read data, models
     *  and other quantities throughout the project.
     * @param NetCDFFile The netcdf file object to read the vector from
     * @param DataName The name of the vector in the netcdf file
     * @param Vec The vector object (std::vector, ublas vector) that will contain the data
     */
    template<class VectorType>
    J3DEXPORT void ReadVec(netCDF::NcFile &NetCDFFile, const std::string &DataName,
        VectorType &Vec)
      {
        // create netcdf variable with the same name as the dimension
        netCDF::NcVar Var = NetCDFFile.getVar(DataName);
        const size_t nvalues = jif3D::cxxport::num_vals(Var);
        //allocate memory in the class variable
        Vec.resize(nvalues);
        //read values from netcdf file

        jif3D::cxxport::get_legacy_ncvar(Var, &Vec[0], nvalues);
      }

    //! Check if a variable exists in a netcdf file
    /*! Check if a variable exists in a netcdf file without reading it.
     *
     * @param NetCDFFile The netcdf file object to check
     * @param Name The name of the variable
     */
    J3DEXPORT inline bool CheckExists(const netCDF::NcFile &NetCDFFile, const std::string &Name)
      {
        bool exists = false;
        try
          {
            netCDF::NcVar Var = NetCDFFile.getVar(Name);
            exists = !Var.isNull();
          } catch (const netCDF::exceptions::NcException &ex)
          {
            exists = false;
          }
        return exists;
      }

    //! Write a vectorial quantity to a netcdf file
    /*! This helper function writes a single vectorial quantity
     * to a netcdf file.
     * @param NetCDFFile The netcdf file object, needs to be writeable
     * @param Name The name of the quantity in the file
     * @param Vec The values to be written
     * @param Dimension The dimension object for the associated dimension, the dimension needs to have identical size to the vector
     * @param unit The units of the values
     */
    template<class VectorType>
    J3DEXPORT void WriteVec(netCDF::NcFile &NetCDFFile, const std::string &Name,
        const VectorType &Vec, const netCDF::NcDim &Dimension, const std::string &unit)
      {
        const size_t nelem = Vec.size();
        netCDF::NcVar PosVar = NetCDFFile.addVar(Name, netCDF::ncDouble, Dimension);
        PosVar.putAtt("units", unit.c_str());

        jif3D::cxxport::put_legacy_ncvar(PosVar, &Vec[0], nelem);
//      PosVar->put(&Vec[0], nelem);
      }

  }

#endif /*NETCDFTOOLS_H_*/

