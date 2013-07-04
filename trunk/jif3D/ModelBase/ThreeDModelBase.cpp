//============================================================================
// Name        : ThreeDModelBase.cpp
// Author      : Max Moorkamp
// Version     :
// Copyright   : 2008, MM
//============================================================================

#include "ThreeDModelBase.h"
#include <cassert>
#include <fstream>
#include <boost/bind.hpp>
#include "../Global/FatalException.h"
#include "NetCDFModelTools.h"
#include "VTKTools.h"
namespace jif3D
  {

    ThreeDModelBase::ThreeDModelBase() :
        MeasPosX(), MeasPosY(), MeasPosZ(), XCellSizesChanged(true), YCellSizesChanged(
            true), ZCellSizesChanged(true), Data(), XCellSizes(), YCellSizes(), ZCellSizes(), GridXCoordinates(), GridYCoordinates(), GridZCoordinates(), XOrigin(
            0.0), YOrigin(0.0), ZOrigin(0.0)
      {
        omp_init_lock(&lck_model_xcoord);
        omp_init_lock(&lck_model_ycoord);
        omp_init_lock(&lck_model_zcoord);

      }

    ThreeDModelBase::~ThreeDModelBase()
      {
        omp_destroy_lock(&lck_model_xcoord);
        omp_destroy_lock(&lck_model_ycoord);
        omp_destroy_lock(&lck_model_zcoord);
      }

    //when we copy a model we always set the cell size change flags to true
    //this costs us recalculation of the grid coordinates, but we do not
    //have to worry whether the grid coordinate values are updated or not
    ThreeDModelBase::ThreeDModelBase(const ThreeDModelBase &source) :
        MeasPosX(source.MeasPosX), MeasPosY(source.MeasPosY), MeasPosZ(source.MeasPosZ), XCellSizesChanged(
            true), YCellSizesChanged(true), ZCellSizesChanged(true), Data(source.Data), XCellSizes(
            source.XCellSizes), YCellSizes(source.YCellSizes), ZCellSizes(
            source.ZCellSizes), GridXCoordinates(), GridYCoordinates(), GridZCoordinates(), XOrigin(
            source.XOrigin), YOrigin(source.YOrigin), ZOrigin(source.ZOrigin)

      {
        //each object needs its own lock for openmp
        //so we do not copy the value from the source object
        //but we reinitialize
        omp_init_lock(&lck_model_xcoord);
        omp_init_lock(&lck_model_ycoord);
        omp_init_lock(&lck_model_zcoord);
      }

    ThreeDModelBase& ThreeDModelBase::operator=(const ThreeDModelBase& source)
      {
        if (this == &source)
          return *this;
        //we have to implement the copy operator to make sure
        //all information is updated appropriately
        //first we copy all the measurement positions
        MeasPosX.resize(source.MeasPosX.size());
        std::copy(source.MeasPosX.begin(), source.MeasPosX.end(), MeasPosX.begin());
        MeasPosY.resize(source.MeasPosY.size());
        std::copy(source.MeasPosY.begin(), source.MeasPosY.end(), MeasPosY.begin());
        MeasPosZ.resize(source.MeasPosZ.size());
        std::copy(source.MeasPosZ.begin(), source.MeasPosZ.end(), MeasPosZ.begin());
        //then we copy the data, i.e. the cells values of the 3D model
        Data.resize(
            boost::extents[source.Data.shape()[0]][source.Data.shape()[1]][source.Data.shape()[2]]);
        Data = source.Data;
        //now we copy the cell sizes for all three directions
        XCellSizes.resize(boost::extents[source.XCellSizes.shape()[0]]);
        std::copy(source.XCellSizes.begin(), source.XCellSizes.end(), XCellSizes.begin());
        YCellSizes.resize(boost::extents[source.YCellSizes.shape()[0]]);
        std::copy(source.YCellSizes.begin(), source.YCellSizes.end(), YCellSizes.begin());
        ZCellSizes.resize(boost::extents[source.ZCellSizes.shape()[0]]);
        std::copy(source.ZCellSizes.begin(), source.ZCellSizes.end(), ZCellSizes.begin());
        //we copy origin of the coordinate system
        XOrigin = source.XOrigin;
        YOrigin = source.YOrigin;
        ZOrigin = source.ZOrigin;
        //we regenerate the coordinate information by pretending that the cell sizes changed
        XCellSizesChanged = true;
        YCellSizesChanged = true;
        ZCellSizesChanged = true;

        //we do not perform any copying of the ompenmp lock
        //they are initialized by the constructor

        return *this;
      }

    /*! This functions assumes that the coordinate of the upper left front corner of the model is
     * is (0,0,0).
     * @param Coordinates This vector will contain the coordinates of the left upper front corner of each cell
     * @param Sizes The size of each cell in m
     * @param ChangeFlag The flag that stores whether this coordinate has been changed
     */
    void ThreeDModelBase::CalcCoordinates(t3DModelDim &Coordinates,
        const t3DModelDim Sizes, bool &ChangeFlag) const
      {
        //create a shorthand for the number of elements
        const size_t nelements = Sizes.size();
        //if the sizes have changed and there is something to calculate
        if (ChangeFlag && nelements > 0)
          {
            //make sure we have enough space for the coordinates
            Coordinates.resize(boost::extents[nelements]);
            //sum up the sizes to get the coordinates
            //in the current setup the first coordinate is always zero
            //this is not ideal and should be changed in the future
            std::partial_sum(Sizes.begin(), Sizes.end() - 1, Coordinates.begin() + 1);
            Coordinates[0] = 0.0;
            ChangeFlag = false;
          }
      }

    boost::array<ThreeDModelBase::t3DModelData::index, 3> ThreeDModelBase::FindAssociatedIndices(
        const double xcoord, const double ycoord, const double zcoord) const
      {
        //for all three directions we use the same approach
        //lower_bound gives an iterator to the field in Coordinates
        //that is just smaller then xcoord
        //the distance between this iterator and the beginning is the index
        const int xindex = std::distance(GetXCoordinates().begin(),
            std::lower_bound(GetXCoordinates().begin(), GetXCoordinates().end(), xcoord));
        const int yindex = std::distance(GetYCoordinates().begin(),
            std::lower_bound(GetYCoordinates().begin(), GetYCoordinates().end(), ycoord));
        const int zindex = std::distance(GetZCoordinates().begin(),
            std::lower_bound(GetZCoordinates().begin(), GetZCoordinates().end(), zcoord));
        //when we return the value we make sure that we cannot go out of bounds
        boost::array<t3DModelData::index, 3> idx =
              {
                { std::max(xindex - 1, 0), std::max(yindex - 1, 0), std::max(zindex - 1,
                    0) } };
        return idx;
      }

    void ThreeDModelBase::SetOrigin(const double x, const double y, const double z)
      {
        //transform the measurement coordinates from old model to new coordinates
        std::transform(MeasPosX.begin(), MeasPosX.end(), MeasPosX.begin(),
            boost::bind(std::plus<double>(), _1, XOrigin - x));
        std::transform(MeasPosY.begin(), MeasPosY.end(), MeasPosY.begin(),
            boost::bind(std::plus<double>(), _1, YOrigin - y));
        std::transform(MeasPosZ.begin(), MeasPosZ.end(), MeasPosZ.begin(),
            boost::bind(std::plus<double>(), _1, ZOrigin - z));
        //copy the information about the new origin
        XOrigin = x;
        YOrigin = y;
        ZOrigin = z;

      }

    void ThreeDModelBase::ReadDataFromNetCDF(const NcFile &NetCDFFile,
        const std::string &DataName, const std::string &UnitsName)
      {
        Read3DModelFromNetCDF(NetCDFFile, DataName, UnitsName, XCellSizes, YCellSizes,
            ZCellSizes, Data);

        //we check that the sizes of the grid cell specifications and the data are matching
        if (XCellSizes.size() != Data.shape()[0] || YCellSizes.size() != Data.shape()[1]
            || ZCellSizes.size() != Data.shape()[2])
          {
            throw jif3D::FatalException("Cell size specification does not match data !");
          }
      }

    void ThreeDModelBase::WriteDataToNetCDF(NcFile &NetCDFFile,
        const std::string &DataName, const std::string &UnitsName) const
      {
        Write3DModelToNetCDF(NetCDFFile, DataName, UnitsName, XCellSizes, YCellSizes,
            ZCellSizes, Data);
      }

    void ThreeDModelBase::WriteVTK(std::string filename,
        const std::string &DataName) const
      {
        Write3DModelToVTK(filename, DataName, GetXCellSizes(), GetYCellSizes(),
            GetZCellSizes(), GetData());
      }

    void ThreeDModelBase::ReadMeasPosNetCDF(const std::string filename)
      {
        jif3D::ReadMeasPosNetCDF(filename, MeasPosX, MeasPosY, MeasPosZ);
      }

    void ThreeDModelBase::ReadMeasPosAscii(const std::string filename)
      {
        //we assume that the measurement positions are simply
        //in the format x,y,z
        std::ifstream infile(filename.c_str());
        double posx, posy, posz;
        while (infile.good())
          {
            infile >> posx >> posy >> posz;
            if (infile.good())
              {
                MeasPosX.push_back(posx);
                MeasPosY.push_back(posy);
                MeasPosZ.push_back(posz);
              }
          }
        assert(MeasPosX.size() == MeasPosY.size());
        assert(MeasPosX.size() == MeasPosZ.size());
      }

  }
