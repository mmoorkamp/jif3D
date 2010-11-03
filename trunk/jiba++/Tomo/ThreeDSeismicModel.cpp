//============================================================================
// Name        : ThreeDSeismicModel.cpp
// Author      : Apr 7, 2009
// Version     :
// Copyright   : 2009, mmoorkamp
//============================================================================


#include "ThreeDSeismicModel.h"
#include "../Global/FatalException.h"
#include <boost/bind.hpp>
#include <cassert>
#include <algorithm>

namespace jiba
  {
    static const std::string SlownessName = "Slowness";
    static const std::string SlownessUnit = "s/m";

    ThreeDSeismicModel::ThreeDSeismicModel() :
      SourcePosX(), SourcePosY(), SourcePosZ(), SourceIndices(),
          ReceiverIndices()
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

    void ThreeDSeismicModel::ReadNetCDF(const std::string filename,
        bool checkgrid)
      {
        //create the netcdf file object
        NcFile DataFile(filename.c_str(), NcFile::ReadOnly);
        //read in the 3D gridded data
        ReadDataFromNetCDF(DataFile, SlownessName, SlownessUnit);
        //check that the grid has equal sizes in all three dimensions
        const double CellSize = GetXCellSizes()[0];
        if (checkgrid)
          {
            if (std::search_n(GetXCellSizes().begin(), GetXCellSizes().end(),
                GetXCellSizes().num_elements(), CellSize)
                == GetXCellSizes().end())
              {
                throw jiba::FatalException(
                    "Non-equal grid spacing in x-direction !");
              }
            if (std::search_n(GetYCellSizes().begin(), GetYCellSizes().end(),
                GetYCellSizes().num_elements(), CellSize)
                == GetYCellSizes().end())
              {
                throw jiba::FatalException(
                    "Non-equal grid spacing in y-direction !");
              }
            if (std::search_n(GetZCellSizes().begin(), GetZCellSizes().end(),
                GetZCellSizes().num_elements(), CellSize)
                == GetZCellSizes().end())
              {
                throw jiba::FatalException(
                    "Non-equal grid spacing in z-direction !");
              }
          }
      }
  }
