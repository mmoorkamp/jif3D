//============================================================================
// Name        : ThreeDSeismicModel.cpp
// Author      : Apr 7, 2009
// Version     :
// Copyright   : 2009, mmoorkamp
//============================================================================

#include "ThreeDSeismicModel.h"
#include "../Global/FatalException.h"
#include <boost/numeric/conversion/cast.hpp>
#include <cassert>
#include <algorithm>
#include <netcdf>

using netCDF::NcFile;
using netCDF::NcVar;
using netCDF::NcDim;

namespace jif3D
  {
    //we use these names when writing the model to a netcdf file
    static const std::string SlownessName = "Slowness";
    static const std::string SlownessUnit = "s/m";

    ThreeDSeismicModel::ThreeDSeismicModel() :
        SourcePosX(), SourcePosY(), SourcePosZ(), SourceIndices(), ReceiverIndices()
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

    ThreeDSeismicModel& ThreeDSeismicModel::operator=(const ThreeDModelBase& source)
      {
        if (this == &source)
          return *this;
        ThreeDModelBase::operator =(source);

        return *this;
      }

    ThreeDSeismicModel& ThreeDSeismicModel::operator=(const ThreeDSeismicModel& source)
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


    boost::array<ThreeDModelBase::t3DModelData::index, 3> ThreeDSeismicModel::FindAssociatedIndices(
        const double xcoord, const double ycoord, const double zcoord) const
      {
        const int xindex = boost::numeric_cast<int>(
            floor((xcoord - GetXCoordinates()[0]) / GetXCellSizes()[0]));
        const int yindex = boost::numeric_cast<int>(
            floor((ycoord - GetYCoordinates()[0]) / GetYCellSizes()[0]));
        const int zindex = boost::numeric_cast<int>(
            floor((zcoord - GetZCoordinates()[0]) / GetZCellSizes()[0]));
        //when we return the value we make sure that we cannot go out of bounds
        boost::array<t3DModelData::index, 3> idx =
          {
            { std::max(xindex, 0), std::max(yindex, 0), std::max(zindex - 1, 0) } };
        return idx;
      }

    void ThreeDSeismicModel::WriteNetCDF(const std::string &filename) const
      {
        NcFile DataFile(filename, NcFile::replace);
        //write the 3D discretized part
        WriteDataToNetCDF(DataFile, SlownessName, SlownessUnit);
      }

    void ThreeDSeismicModel::ReadNetCDF(const std::string &filename)
      {
        ReadNetCDF(filename, false);
      }

    void ThreeDSeismicModel::ReadNetCDF(const std::string &filename, bool checkgrid)
      {
        //create the netcdf file object
        NcFile DataFile(filename, NcFile::read);
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
                GetXCellSizes().size(), CellSize) == GetXCellSizes().end())
              {
                throw jif3D::FatalException("Non-equal grid spacing in x-direction !",
                __FILE__, __LINE__);
              }
            //then for y
            if (std::search_n(GetYCellSizes().begin(), GetYCellSizes().end(),
                GetYCellSizes().size(), CellSize) == GetYCellSizes().end())
              {
                throw jif3D::FatalException("Non-equal grid spacing in y-direction !",
                __FILE__, __LINE__);
              }
            //finally for z, in each cases the cell size we search for is the same
            if (std::search_n(GetZCellSizes().begin(), GetZCellSizes().end(),
                GetZCellSizes().size(), CellSize) == GetZCellSizes().end())
              {
                throw jif3D::FatalException("Non-equal grid spacing in z-direction !",
                __FILE__, __LINE__);
              }
          }
      }
  }
