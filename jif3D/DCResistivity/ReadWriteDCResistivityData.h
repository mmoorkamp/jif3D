//============================================================================
// Name        : ReadWriteDCResistivityData.h
// Author      : May 8, 2014
// Version     :
// Copyright   : 2014, zhanjie and mmoorkamp
//============================================================================

#ifndef READWRITEDCRESISTIVITYDATA_H_
#define READWRITEDCRESISTIVITYDATA_H_

#include <fstream>
#include "../Global/VecMat.h"
#include "ThreeDDCResistivityModel.h"


namespace jif3D
  {
    /** \addtogroup dcresistivity classes and functions */
    /* @{ */
    //!Save a collection of apparent resistivity data and associated positions to a netcdf file
    /*! Saves apparent resistivity data together with the source receiver configuration to a netcdf file.
     * @param filename The name of the file to store the data
     * @param Data The apparent resistivity data in ohm.m, the length of this vector matches the number of
     *             source-receiver combinations in the Model object
     * @param The error of the apparent resistivity
     * @param Model The model object containing information about the source and receiver setup
     */
    void SaveApparentResistivity(const std::string &filename, const jif3D::rvec &Data,
        const jif3D::rvec &Error, const jif3D::ThreeDDCResistivityModel &Model);
    //!Read a collection of apparent resistivity and associated positions from a netcdf file
    /*! Read apparent resistivity together with the source receiver configuration from a netcdf file.
     * @param filename The name of the file to store the data
     * @param Data The apparent resistivity in ohm.m, the length of this vector matches the number of
     *         source-receiver combinations in the Model object
     * @param The error of the apparent resistivity
     * @param Model The model object containing information about the source and receiver setup
     */
    void ReadApparentResistivity(const std::string &filename, jif3D::rvec &Data,
        jif3D::rvec &Error, jif3D::ThreeDDCResistivityModel &Model);

  /* @} */
  }
#endif /* READWRITETOMOGRAPHYDATA_H_ */
