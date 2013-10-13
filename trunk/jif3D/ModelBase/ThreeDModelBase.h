//============================================================================
// Name        : ThreeDModelBase.h
// Author      : Max Moorkamp
// Version     :
// Copyright   : 2008, MM
//============================================================================

#ifndef THREEDMODELBASE_H_
#define THREEDMODELBASE_H_

#include <netcdfcpp.h>
#ifdef HAVEOPENMP
#include <omp.h>
#endif
#include <string>
#include <boost/multi_array.hpp>
#include <boost/serialization/serialization.hpp>
#include <boost/serialization/vector.hpp>

/*! \file ThreeDModelBase.h
 * Contains the base class for all 3D models.
 */
namespace jif3D
  {
    /** \addtogroup modelbase Basic classes and routines for 3D models
     * This module contains classes that are associated with 3D
     * rectilinear gridded models in general without restriction to a certain method.
     * This unified approach allows a transparent
     * handling of model information within all codes. Furthermore it is easy
     * to convert a model for one method into a model for another method, e.g.
     * to realize a fixed relationship between model parameters.
     */
    /* @{ */

    //! The basic storage class for three-dimensional models
    /*! This class provides the storage and general functionality
     * for types all of three-dimensional models. Any other class
     * that needs to store a three-dimensional model should be derived from it.
     */
    class ThreeDModelBase
      {
    public:
      //! A shorthand for the type of the stored model data for 3D models
      typedef boost::multi_array<double, 3> t3DModelData;
      //! A shorthand for the dimensions, i.e. cell-sizes for the 3D models
      typedef boost::multi_array<double, 1> t3DModelDim;
      //! The type of the measurement position vector, this is a std::vector because we want to easily append elements
      typedef std::vector<double> tMeasPosVec;
    private:
      //! the x-coordinates of the measurement points
      tMeasPosVec MeasPosX;
      //! the y-coordinates of the measurement points
      tMeasPosVec MeasPosY;
      //! the z-coordinates of the measurement points
      tMeasPosVec MeasPosZ;
      //! Have the cell sizes for the x-coordinate changed
      /*! This is basically a flag for caching, we could recalculate the Grid coordinates
       * every time from the cell sizes. This is slow and the coordinates are often needed.
       * Therefore we store whether the sizes have changed and only recalculate then, otherwise
       * we take the values from GridXCoordinates. This should also be possible for const
       * objects, even though the object is not bitwise constant any more. Thus it is mutable.
       */
      mutable bool XCellSizesChanged;
      //! Have the cell sizes for the y-coordinate changed
      mutable bool YCellSizesChanged;
      //! Have the cell sizes for the z-coordinate changed
      mutable bool ZCellSizesChanged;
#ifdef HAVEOPENMP
      //! a locking variable to allow concurrent coordinate calculation calls
      mutable omp_lock_t lck_model_coord;
#endif
      //! The object containing the actual value, e.g. conductivity, velocity
      t3DModelData Data;
      //! The size of the cells in x-direction
      t3DModelDim XCellSizes;
      //! The size of the cells in y-direction
      t3DModelDim YCellSizes;
      //! The size of the cells in z-direction
      t3DModelDim ZCellSizes;
      //! The x-coordinate of the upper left front corner
      /*! See the explanation for XCellSizesChanged why this is mutable
       */
      mutable t3DModelDim GridXCoordinates;
      //! The y-coordinate of the upper left front corner
      mutable t3DModelDim GridYCoordinates;
      //! The z-coordinate of the upper left front corner
      mutable t3DModelDim GridZCoordinates;
      //! Calculate the coordinates of the model cells from the sizes of each cell, this is a helper function used by ThreeDModelBase
      void CalcCoordinates(t3DModelDim &Coordinates, const t3DModelDim Sizes,
          bool &ChangeFlag) const;
      friend class boost::serialization::access;
      //! Provide serialization to be able to store objects
      template<class Archive>
      void save(Archive & ar, const unsigned int version) const
        {
          //we can directly serialize std::vector, so MeasPos are simple
          ar & MeasPosX;
          ar & MeasPosY;
          ar & MeasPosZ;
          //multi-array does not have serialization support
          //so we have to store the shape of the arrays
          ar & Data.shape()[0];
          ar & Data.shape()[1];
          ar & Data.shape()[2];
          //then serialize the raw data
          ar & boost::serialization::make_array(Data.origin(), Data.num_elements());
          ar
              & boost::serialization::make_array(XCellSizes.origin(),
                  XCellSizes.num_elements());
          ar
              & boost::serialization::make_array(YCellSizes.origin(),
                  YCellSizes.num_elements());
          ar
              & boost::serialization::make_array(ZCellSizes.origin(),
                  ZCellSizes.num_elements());
          //finally we need the origin of the model
          ar & XOrigin;
          ar & YOrigin;
          ar & ZOrigin;
        }
      //! Provide serialization to be able to load objects
      template<class Archive>
      void load(Archive & ar, const unsigned int version)
        {
          //we can directly serialize std::vector, so MeasPos are simple
          ar & MeasPosX;
          ar & MeasPosY;
          ar & MeasPosZ;
          //multi-array does not have serialization support
          //so we have to load the shape of the arrays
          size_t nx, ny, nz;
          ar & nx;
          ar & ny;
          ar & nz;
          Data.resize(boost::extents[nx][ny][nz]);
          XCellSizes.resize(boost::extents[nx]);
          YCellSizes.resize(boost::extents[ny]);
          ZCellSizes.resize(boost::extents[nz]);
          ar & boost::serialization::make_array(Data.origin(), Data.num_elements());
          ar
              & boost::serialization::make_array(XCellSizes.origin(),
                  XCellSizes.num_elements());
          ar
              & boost::serialization::make_array(YCellSizes.origin(),
                  YCellSizes.num_elements());
          ar
              & boost::serialization::make_array(ZCellSizes.origin(),
                  ZCellSizes.num_elements());
          ar & XOrigin;
          ar & YOrigin;
          ar & ZOrigin;
        }
      BOOST_SERIALIZATION_SPLIT_MEMBER()protected:
      //! The origin of the coordinate system in x-direction in m
      double XOrigin;
      //! The origin of the coordinate system in y-direction in m
      double YOrigin;
      //! The origin of the coordinate system in z-direction in m
      double ZOrigin;
      //! read-write access to the cell size in x-direction in m
      t3DModelDim &SetXCellSizes()
        {
          XCellSizesChanged = true;
          return XCellSizes;
        }
      //! read-write access to the cell size in y-direction in m
      t3DModelDim &SetYCellSizes()
        {
          YCellSizesChanged = true;
          return YCellSizes;
        }
      //! read-write access to the cell size in z-direction in m
      t3DModelDim &SetZCellSizes()
        {
          ZCellSizesChanged = true;
          return ZCellSizes;
        }
      //! Read data and associated cell sizes from a netcdf file
      void ReadDataFromNetCDF(const NcFile &NetCDFFile, const std::string &DataName,
          const std::string &UnitsName);
      //! Write data and associated cell sizes to a netcdf file
      void WriteDataToNetCDF(NcFile &NetCDFFile, const std::string &DataName,
          const std::string &UnitsName) const;
      //! Write the data and cell sizes to a VTK file for plotting in Paraview or Visit etc.
      void WriteVTK(std::string filename, const std::string &DataName) const;
    public:
      //! return a reference to the data so it can be modified
      t3DModelData &SetData()
        {
          return Data;
        }
      //! return read only access to the stored data
      const t3DModelData &GetData() const
        {
          return Data;
        }
      //! Add a measurement point to the model
      void AddMeasurementPoint(const double xcoord, const double ycoord,
          const double zcoord)
        {
          MeasPosX.push_back(xcoord - XOrigin);
          MeasPosY.push_back(ycoord - YOrigin);
          MeasPosZ.push_back(zcoord - ZOrigin);
        }
      //! remove all information about measurement points
      void ClearMeasurementPoints()
        {
          MeasPosX.clear();
          MeasPosY.clear();
          MeasPosZ.clear();
        }
      //! Return the x-coordinates (Northing) of all measurement points read-only
      /*! This function provides read-only access to the x-coordinates
       * of the measurement points. The only way to modify the position of
       * the measurements is to delete them with ClearMeasurementPoints and
       * add new ones with AddMeasurementPoint. This ensures that we have all
       * three coordinate values for all points.
       * @return A vector with the x-coordinates of all measurement points in m
       */
      const tMeasPosVec &GetMeasPosX() const
        {
          return MeasPosX;
        }
      //! Return the y-coordinates (Easting)of all measurement points read-only
      const tMeasPosVec &GetMeasPosY() const
        {
          return MeasPosY;
        }
      //! Return the z-coordinates (Depth) of all measurement points read-only
      const tMeasPosVec &GetMeasPosZ() const
        {
          return MeasPosZ;
        }
      //! Set the origin of the coordinate system
      virtual void SetOrigin(const double x, const double y, const double z);
      //! Set the size of the mesh and the coordinate axes
      void SetMeshSize(const size_t nx, const size_t ny, const size_t nz)
        {
          XCellSizes.resize(boost::extents[nx]);
          YCellSizes.resize(boost::extents[ny]);
          ZCellSizes.resize(boost::extents[nz]);
          Data.resize(boost::extents[nx][ny][nz]);
        }
      //! From the three spatial indices, calculate the offset in memory
      int IndexToOffset(int xi, int yi, int zi) const
        {
          return Data.shape()[2] * (xi * Data.shape()[1] + yi) + zi;
        }
      //! Form a memory offset, calculate the associated spatial indices
      void OffsetToIndex(int offset, int &xi, int &yi, int &zi) const
        {
          zi = offset % Data.shape()[2];
          xi = (offset - zi) / Data.shape()[2];
          yi = xi % Data.shape()[1];
          xi = (xi - yi) / Data.shape()[1];
        }
      //! Return the size of the gridded domain in x, y and z-direction, respectively
      /*! This function solves the problem that sometimes we want to know how
       * big the gridded domain is for a general model (seismic, MT or gravity)
       * but Data is a private member. This function just returns
       * The result of Data.shape()
       * @return The size of the gridded domain in x, y and z-direction as it would be returned by a call to Data.shape()
       */
      const boost::multi_array_types::size_type* GetModelShape() const
        {
          return Data.shape();
        }
      //! Return the total number of cells in the gridded domain
      size_t GetNModelElements() const
        {
          return Data.num_elements();
        }
      //! read-only access to the cell size in x-direction in m
      const t3DModelDim &GetXCellSizes() const
        {
          return XCellSizes;
        }

      //! read-only access to the cell size in y-direction in m
      const t3DModelDim &GetYCellSizes() const
        {
          return YCellSizes;
        }

      //! read-only access to the cell size in z-direction in m
      const t3DModelDim &GetZCellSizes() const
        {
          return ZCellSizes;
        }

      //!Get the x (north) coordinates of the cells, might perform calculations and write operation but is now thread-safe
      const t3DModelDim &GetXCoordinates() const
        {
          CalcCoordinates(GridXCoordinates, XCellSizes, XCellSizesChanged);
          return GridXCoordinates;
        }
      //!Get the y (east) coordinates of the cells, might perform calculations and write operation but is now thread-safe
      const t3DModelDim &GetYCoordinates() const
        {
          CalcCoordinates(GridYCoordinates, YCellSizes, YCellSizesChanged);
          return GridYCoordinates;
        }
      //!Get the z (depth) coordinates of the cells, might perform calculations and write operation but is now thread-safe
      const t3DModelDim &GetZCoordinates() const
        {
          CalcCoordinates(GridZCoordinates, ZCellSizes, ZCellSizesChanged);
          return GridZCoordinates;
        }
      //! The derived model classes also manage the synthetic data configuration, for parallel calculations we can signal how many chunks can be calculated independently (e.g. how many frequencies)
      virtual size_t GetNIndependentChunks() const
        {
          return 1;
        }
      //! Given three coordinates in m, find the indices of the model cell that correponds to these coordinates
      virtual boost::array<ThreeDModelBase::t3DModelData::index, 3>
      FindAssociatedIndices(const double xcoord, const double ycoord,
          const double zcoord) const;
      //! Read the Measurement positions from a netcdf file
      void ReadMeasPosNetCDF(const std::string filename);
      //! Read the Measurement positions from an ascii file
      void ReadMeasPosAscii(const std::string filename);
      friend class ModelRefiner;
      ThreeDModelBase();
      //! The copy operator copies all independent variables of ThreeDModelBase
      ThreeDModelBase& operator=(const ThreeDModelBase& source);
      //! We need to define a copy constructor to deal with openmp locks
      ThreeDModelBase(const ThreeDModelBase &source);
      virtual ~ThreeDModelBase();
    };
/* @} */
}

#endif /*THREEDMODELBASE_H_*/
