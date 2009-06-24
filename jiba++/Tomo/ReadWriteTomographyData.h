//============================================================================
// Name        : ReadWriteTomographyData.h
// Author      : May 5, 2009
// Version     :
// Copyright   : 2009, mmoorkamp
//============================================================================


#ifndef READWRITETOMOGRAPHYDATA_H_
#define READWRITETOMOGRAPHYDATA_H_

#include "../Global/VecMat.h"
#include "ThreeDSeismicModel.h"
#include "modeling_seismic.h"

namespace jiba
  {
    /** \addtogroup tomo Seismic tomography classes and functions */
    /* @{ */
    //!Save a collection of travel times and associated positions to a netcdf file
    void SaveTraveltimes(const std::string &filename, const jiba::rvec &Data,
        const jiba::ThreeDSeismicModel &Model);
    //!Read a collection of travel times and associated positions from a netcdf file
    void ReadTraveltimes(const std::string &filename, jiba::rvec &Data,
        jiba::ThreeDSeismicModel &Model);
    //! Plot the travel time field, i.e. the travel time at each grid cell in a netcdf file
    void PlotTimeField(const std::string &filename, const float *Times,
        const double gridspacing, const size_t nx, const size_t ny,
        const size_t nz);
    //! Plot the raypaths for a 3D forward calculation in a .vtk file
    void PlotRaypath(const std::string &filename, jiba::RP_STRUCT *raypath,
        const size_t nmeas, const double gridspacing, const size_t nairlayers);
  /* @} */
  }
#endif /* READWRITETOMOGRAPHYDATA_H_ */
