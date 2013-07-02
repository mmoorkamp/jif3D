//============================================================================
// Name        : NetCDFTools.h
// Author      : 1 Jul 2013
// Version     : 
// Copyright   : 2013, mm489
//============================================================================

#ifndef NETCDFTOOLS_H_
#define NETCDFTOOLS_H_

#include <netcdfcpp.h>

namespace jif3D {

  //! Read a vector from a netcdf file
  /*! This helper function reads a single vectorial quantity from
   *  an already opened netcdf file. It is used to read data, models
   *  and other quantities throughout the project.
   * @param NetCDFFile The netcdf file object to read the vector from
   * @param DataName The name of the vector in the netcdf file
   * @param Vec The vector object (std::vector, ublas vector) that will contain the data
   */
  template<class VectorType>
  void ReadVec(NcFile &NetCDFFile, const std::string &DataName,
      VectorType &Vec)
    {
    // create netcdf variable with the same name as the dimension
      NcVar *Var = NetCDFFile.get_var(DataName.c_str());
      const size_t nvalues = Var->num_vals();
      //allocate memory in the class variable
      Vec.resize(nvalues);
      //read values from netcdf file
      Var->get(&Vec[0], nvalues);
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
  void WriteVec(NcFile &NetCDFFile, const std::string &Name,
      const VectorType &Vec, NcDim *Dimension, const std::string unit)
    {
      const size_t nelem = Vec.size();
      NcVar *PosVar = NetCDFFile.add_var(Name.c_str(), ncDouble,
          Dimension);
      PosVar->add_att("units", unit.c_str());
      PosVar->put(&Vec[0], nelem);
    }

}

#endif /*NETCDFTOOLS_H_*/

