//============================================================================
// Name        : ThreeDSeismicModel.cpp
// Author      : Apr 7, 2009
// Version     :
// Copyright   : 2009, mmoorkamp
//============================================================================


#include "ThreeDSeismicModel.h"
#include <boost/bind.hpp>
#include <cassert>

namespace jiba
  {
    static const std::string SlownessName = "Slowness";
    static const std::string SlownessUnit = "s/km";

    ThreeDSeismicModel::ThreeDSeismicModel()
      {
      }

    ThreeDSeismicModel::~ThreeDSeismicModel()
      {
      }

    void ThreeDSeismicModel::SetOrigin(const double x, const double y,
        const double z)
      {
        //transform the source coordinates from old model to real coordinates
        std::transform(SourcePosX.begin(), SourcePosX.end(),
            SourcePosX.begin(), boost::bind(std::plus<double>(), _1, XOrigin
                - x));
        std::transform(SourcePosY.begin(), SourcePosY.end(),
            SourcePosY.begin(), boost::bind(std::plus<double>(), _1, YOrigin
                - y));
        std::transform(SourcePosZ.begin(), SourcePosZ.end(),
            SourcePosZ.begin(), boost::bind(std::plus<double>(), _1, ZOrigin
                - z));
        //we have to call the base implementation in the end because
        //it changes the measurement positions and the Origin
        ThreeDModelBase::SetOrigin(x, y, z);
      }

    void ThreeDSeismicModel::WriteNetCDF(const std::string filename) const
      {

        NcFile DataFile(filename.c_str(), NcFile::Replace);
        //first write the 3D discretized part
        WriteDataToNetCDF(DataFile, SlownessName, SlownessUnit);

      }

    void ThreeDSeismicModel::ReadNetCDF(const std::string filename)
      {
        //create the netcdf file object
        NcFile DataFile(filename.c_str(), NcFile::ReadOnly);
        //read in the 3D gridded data
        ReadDataFromNetCDF(DataFile, SlownessName, SlownessUnit);
      }
  }
