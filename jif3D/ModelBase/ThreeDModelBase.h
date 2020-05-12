//============================================================================
// Name        : ThreeDModelBase.h
// Author      : Max Moorkamp
// Version     :
// Copyright   : 2008, MM
//============================================================================

#ifndef THREEDMODELBASE_H_
#define THREEDMODELBASE_H_

#include "../Global/Serialization.h"
#include "../Global/FatalException.h"
#include "../Global/Jif3DGlobal.h"
#include <string>
#include <netcdf>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wmisleading-indentation"
#pragma GCC diagnostic ignored "-Wunused-local-typedefs"
#include <boost/multi_array.hpp>
#pragma GCC diagnostic pop

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
    class J3DEXPORT ThreeDModelBase
      {
    public:
      friend class MatOpRegularization;
      //! A shorthand for the type of the stored model data for 3D models
      typedef boost::multi_array<double, 3> t3DModelData;
      //! A shorthand for the dimensions, i.e. cell-sizes for the 3D models
      typedef std::vector<double> t3DModelDim;
    private:

      //! The object containing the actual value, e.g. conductivity, velocity
      t3DModelData Data;
      //! The size of the cells in x-direction
      t3DModelDim XCellSizes;
      //! The size of the cells in y-direction
      t3DModelDim YCellSizes;
      //! The size of the cells in z-direction
      t3DModelDim ZCellSizes;
      //! The coordinate of the grid cells in x-direction, contains both corners and therefore has nx+1 elements
      t3DModelDim GridXCoordinates;
      //! The coordinate of the grid cells in y-direction, contains both corners and therefore has ny+1 elements
      t3DModelDim GridYCoordinates;
      //! The coordinate of the grid cells in y-direction, contains both corners and therefore has nz+1 elements
      t3DModelDim GridZCoordinates;
      //! Calculate the coordinates of the model cells from the sizes of each cell, this is a helper function used by ThreeDModelBase
      void CalcSizes(const t3DModelDim &Coordinates, t3DModelDim &Sizes);
      void CalcCoords(t3DModelDim &Coordinates, const t3DModelDim &Sizes);
    protected:
      //! write access to the cell size in x-direction in m
      void SetXCellSizes(const t3DModelDim &XCS)
        {
          XCellSizes = XCS;
          CalcCoords(GridXCoordinates, XCellSizes);
        }
      //! write access to the cell size in y-direction in m
      void SetYCellSizes(const t3DModelDim &YCS)
        {
          YCellSizes = YCS;
          CalcCoords(GridYCoordinates, YCellSizes);
        }
      //! write access to the cell size in z-direction in m
      void SetZCellSizes(const t3DModelDim &ZCS)
        {
          ZCellSizes = ZCS;
          CalcCoords(GridZCoordinates, ZCellSizes);
        }
      //! write access to the cell size in x-direction in m
      void SetXCoordinates(const t3DModelDim &XCS)
        {
          GridXCoordinates = XCS;
          CalcSizes(GridXCoordinates, XCellSizes);
        }
      //! write access to the cell size in y-direction in m
      void SetYCoordinates(const t3DModelDim &YCS)
        {
          GridYCoordinates = YCS;
          CalcSizes(GridYCoordinates, YCellSizes);
        }
      //! write access to the cell size in z-direction in m
      void SetZCoordinates(const t3DModelDim &ZCS)
        {
          GridZCoordinates = ZCS;
          CalcSizes(GridZCoordinates, ZCellSizes);
        }
      //! Read data and associated cell sizes from a netcdf file
      void ReadDataFromNetCDF(const netCDF::NcFile &NetCDFFile,
          const std::string &DataName, const std::string &UnitsName);
      //! Write data and associated cell sizes to a netcdf file
      void WriteDataToNetCDF(netCDF::NcFile &NetCDFFile, const std::string &DataName,
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
      //! Set the origin of the coordinate system
      virtual void SetOrigin(const double x, const double y, const double z);
      //! Get the Model origin for the x-coordinate
      double GetXOrigin() const
        {
          return GridXCoordinates.empty() ? 0.0 : GridXCoordinates.at(0);
        }
      //! Get the Model origin for the y-coordinate
      double GetYOrigin() const
        {
          return GridYCoordinates.empty() ? 0.0 : GridYCoordinates.at(0);
        }
      //! Get the Model origin for the z-coordinate
      double GetZOrigin() const
        {
          return GridZCoordinates.empty() ? 0.0 : GridZCoordinates.at(0);
        }
      //! Set the size of the mesh and the coordinate axes
      void SetMeshSize(const size_t nx, const size_t ny, const size_t nz)
        {
          XCellSizes.resize(nx);
          YCellSizes.resize(ny);
          ZCellSizes.resize(nz);
          GridXCoordinates.resize(nx + 1);
          GridYCoordinates.resize(ny + 1);
          GridZCoordinates.resize(nz + 1);
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
          return GridXCoordinates;
        }
      //!Get the y (east) coordinates of the cells, might perform calculations and write operation but is now thread-safe
      const t3DModelDim &GetYCoordinates() const
        {
          return GridYCoordinates;
        }
      //!Get the z (depth) coordinates of the cells, might perform calculations and write operation but is now thread-safe
      const t3DModelDim &GetZCoordinates() const
        {
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
      virtual void WriteNetCDF(const std::string &filename) const
        {
          throw jif3D::FatalException("WriteNetCDF not implemented in ThreeDModelBase !");
        }
      //! Read all model information from a netcdf file
      virtual void ReadNetCDF(const std::string &filename)
        {
          throw jif3D::FatalException("ReadNetCDF not implemented in ThreeDModelBase !");
        }
      //! Write model in the format x, y, z, value, where x/y/z are the centres of the cells
      void WriteXYZ(const std::string &filename) const;
      //! Write model in the format used by UBC and Geoscience Analyst
      void WriteUBC(const std::string &filename) const;
      //! Provide serialization to be able to store objects
      template<class Archive>
      void save(Archive & ar, const unsigned int version) const
        {
          //multi-array does not have serialization support
          //so we have to store the shape of the arrays
          ar & Data.shape()[0];
          ar & Data.shape()[1];
          ar & Data.shape()[2];
          //then serialize the raw data
          ar & make_array(Data.origin(), Data.num_elements());
          ar & GridXCoordinates;
          ar & GridYCoordinates;
          ar & GridZCoordinates;

        }
      //! Provide serialization to be able to load objects
      template<class Archive>
      void load(Archive & ar, const unsigned int version)
        {
          //multi-array does not have serialization support
          //so we have to load the shape of the arrays
          size_t nx = 0, ny = 0, nz = 0;
          ar & nx;
          ar & ny;
          ar & nz;
          Data.resize(boost::extents[nx][ny][nz]);
          ar & make_array(Data.origin(), Data.num_elements());
          ar & GridXCoordinates;
          ar & GridYCoordinates;
          ar & GridZCoordinates;
          CalcSizes(GridXCoordinates, XCellSizes);
          CalcSizes(GridYCoordinates, YCellSizes);
          CalcSizes(GridZCoordinates, ZCellSizes);
        }
#ifdef HAVEHPX
      HPX_SERIALIZATION_SPLIT_MEMBER()
#else
      BOOST_SERIALIZATION_SPLIT_MEMBER()
#endif
      friend class ModelRefiner;
      ThreeDModelBase();
      ThreeDModelBase(const ThreeDModelBase& source);
      bool operator ==(const ThreeDModelBase &b) const;
      ThreeDModelBase& operator=(const ThreeDModelBase& source);
      virtual ~ThreeDModelBase();
      };
  /* @} */
  }

#endif /*THREEDMODELBASE_H_*/
