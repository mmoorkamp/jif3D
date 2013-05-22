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

namespace jif3D
  {
    //we use these names when writing the model to a netcdf file
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

    ThreeDSeismicModel::ThreeDSeismicModel(const ThreeDSeismicModel &source) :
      ThreeDModelBase(source), SourcePosX(source.SourcePosX), SourcePosY(
          source.SourcePosY), SourcePosZ(source.SourcePosZ), SourceIndices(
          source.SourceIndices), ReceiverIndices(source.ReceiverIndices)
      {

      }

    ThreeDSeismicModel& ThreeDSeismicModel::operator=(
        const ThreeDSeismicModel& source)
      {
        if (this == &source)
          return *this;
        ThreeDModelBase::operator =(source);
        SourcePosX = source.SourcePosX;
        SourcePosY = source.SourcePosY;
        SourcePosZ = source.SourcePosZ;
        SourceIndices = source.SourceIndices;
        ReceiverIndices = source.ReceiverIndices;
        return *this;
      }

    void ThreeDSeismicModel::SetOrigin(const double x, const double y,
        const double z)
      {
        //transform the source coordinates from old model to real coordinates
        //the coordinates of the receivers are changed by the implementation
        //in the base class that we call below
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
        //write the 3D discretized part
        WriteDataToNetCDF(DataFile, SlownessName, SlownessUnit);

      }

    void ThreeDSeismicModel::ReadNetCDF(const std::string filename,
        bool checkgrid)
      {
        //create the netcdf file object
        NcFile DataFile(filename.c_str(), NcFile::ReadOnly);
        //read in the 3D gridded data
        ReadDataFromNetCDF(DataFile, SlownessName, SlownessUnit);

        const double CellSize = GetXCellSizes()[0];
        //we can check that the grid has equal sizes in all three dimensions
        //in some cases we use a seismic grid that does not conform to
        //the requirements of the forward code, e.g. as an inversion grid
        //that is then refined for the forward calculation, that is why
        //we can turn of checking through the parameter checkgrid
        if (checkgrid)
          {
            //all grid cells should have the same size, so we search for
            // n = number of cells occurences of CellSize, if the search
            //fails search_n returns the end iterator and we throw an exception
            //we do this for all three spatial directions
            //first for x
            if (std::search_n(GetXCellSizes().begin(), GetXCellSizes().end(),
                GetXCellSizes().num_elements(), CellSize)
                == GetXCellSizes().end())
              {
                throw jif3D::FatalException(
                    "Non-equal grid spacing in x-direction !");
              }
            //then for y
            if (std::search_n(GetYCellSizes().begin(), GetYCellSizes().end(),
                GetYCellSizes().num_elements(), CellSize)
                == GetYCellSizes().end())
              {
                throw jif3D::FatalException(
                    "Non-equal grid spacing in y-direction !");
              }
            //finally for z, in each cases the cell size we search for is the same
            if (std::search_n(GetZCellSizes().begin(), GetZCellSizes().end(),
                GetZCellSizes().num_elements(), CellSize)
                == GetZCellSizes().end())
              {
                throw jif3D::FatalException(
                    "Non-equal grid spacing in z-direction !");
              }
          }
      }
  }
