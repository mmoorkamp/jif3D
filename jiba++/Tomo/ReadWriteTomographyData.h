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
#include "inv3d.h"

namespace jiba
  {
    //!Save a collection of scalar gravity measurements and associated positions to a netcdf file
    void SaveTraveltimes(const std::string &filename, const jiba::rvec &Data,
        const jiba::ThreeDSeismicModel &Model);
    void ReadTraveltimes(const std::string &filename, jiba::rvec &Data,
        jiba::ThreeDSeismicModel &Model);
    void PlotTimeField(const std::string &filename, const float *Times,
        const double gridspacing, const size_t nx, const size_t ny,
        const size_t nz);
    void PlotRaypath(const std::string &filename, jiba::RP_STRUCT *raypath,
        const size_t nmeas, const double gridspacing, const size_t nairlayers);
  }
#endif /* READWRITETOMOGRAPHYDATA_H_ */
