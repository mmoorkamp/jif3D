//============================================================================
// Name        : ThreeDModelBase.h
// Author      : Max Moorkamp
// Version     :
// Copyright   : 2008, MM
//============================================================================


#ifndef THREEDMODELBASE_H_
#define THREEDMODELBASE_H_
#include <boost/multi_array.hpp>
#include <netcdfcpp.h>
#include <string>
#include <functional>
#include <boost/bind.hpp>
namespace jiba
  {
    //! The basic storage classe for three-dimensional models
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
  private:
      //! If a derived class needs to perform some action when one the cell size changes, this can be done here
      virtual void SetCellSizesAction(t3DModelDim &sizes)
        {
        }
      //! Have the cell sizes for the x-coordinate changed
      bool XCellSizesChanged;
      //! Have the cell sizes for the y-coordinate changed
      bool YCellSizesChanged;
      //! Have the cell sizes for the z-coordinate changed
      bool ZCellSizesChanged;
      //! The object containing the actual value, e.g. conductivity, velocity
      t3DModelData Data;
      //! The size of the cells in x-direction
      t3DModelDim XCellSizes;
      //! The size of the cells in y-direction
      t3DModelDim YCellSizes;
      //! The size of the cells in z-direction
      t3DModelDim ZCellSizes;
      //! The x-coordinate of the upper left front corner
      t3DModelDim GridXCoordinates;
      //! The x-coordinate of the upper left front corner
      t3DModelDim GridYCoordinates;
      //! The x-coordinate of the upper left front corner
      t3DModelDim GridZCoordinates;
      //! Read a single CellSize variable from a netcdf file
      /*! This is a helper function that reads the cell length for a single
       * dimension from the file.
       */
      size_t ReadSizesFromNetCDF(const NcFile &NetCDFFile,
          const std::string &SizeName, t3DModelDim &CellSize,
          t3DModelDim CellCoordinates);
      //! Write the length of the model cell along a single dimension to the file
      NcDim *WriteSizesToNetCDF(NcFile &NetCDFFile,
          const std::string &SizeName, const t3DModelDim &CellSize) const;
      //! Calculate the coordinates of the model cells from the sizes of each cell
      /*! This functions assumes that the coordinate of the upper left front corner of the model is
       * is (0,0,0).
       * @param Coordinates This vector will contain the coordinates of the left upper front corner of each cell
       * @param Sizes The size of each cell in m
       * @param ChangeFlag The flag that stores whether this coordinate has been changed
       */
      void CalcCoordinates(t3DModelDim &Coordinates, const t3DModelDim Sizes,
          bool &ChangeFlag)
        {
          Coordinates.resize(boost::extents[Sizes.size()]);
          if (Sizes.size() > 0)
            {
              std::partial_sum(Sizes.begin(), Sizes.end(), Coordinates.begin());
              std::rotate(Coordinates.begin(), Coordinates.end()-1,
                  Coordinates.end());
              Coordinates[0] = 0.0;
            }
          ChangeFlag = false;
        }
  protected:
      //The access to the data is protected, this requires the derived class
      // to provide an access function and allows to add initial checks and
      // meaningful names
      //! return read only access to the stored data
      const t3DModelData &GetData() const
        {
          return Data;
        }
      //! return a reference to the data so it can be modified
      t3DModelData &SetData()
        {
          return Data;
        }
      //!Get the x (north) coordinates of the cells
      const t3DModelDim &GetXCoordinates()
        {
          if (XCellSizesChanged && XCellSizes.size() > 0)
            {
              CalcCoordinates(GridXCoordinates,XCellSizes,XCellSizesChanged);
            }
          return GridXCoordinates;
        }
      //!Get the y (east) coordinates of the cells
      const t3DModelDim &GetYCoordinates()
        {
          if (YCellSizesChanged && YCellSizes.size() > 0)
            {
              CalcCoordinates(GridYCoordinates,YCellSizes,YCellSizesChanged);
            }
          return GridYCoordinates;
        }
      //!Get the z (depth) coordinates of the cells
      const t3DModelDim &GetZCoordinates()
        {
          if (ZCellSizesChanged && ZCellSizes.size() > 0)
            {
              CalcCoordinates(GridZCoordinates,ZCellSizes,ZCellSizesChanged);
            }
          return GridZCoordinates;
        }
      //! Read data and associated cell sizes from a netcdf file
      void ReadDataFromNetCDF(const NcFile &NetCDFFile,
          const std::string &DataName, const std::string &UnitsName);
      //! Write data and associated cell sizes to a netcdf file
      void WriteDataToNetCDF(NcFile &NetCDFFile, const std::string &DataName,
          const std::string &UnitsName) const;
      //! Write the data and cell sizes to a VTK file for plotting in Paraview or Visit etc.
      void WriteVTK(std::string filename, const std::string &DataName);
  public:
      //! read-only access to the cell size in x-direction in m
      const t3DModelDim &GetXCellSizes() const
        {
          return XCellSizes;
        }
      //! read-write access to the cell size in x-direction in m
      t3DModelDim &SetXCellSizes()
        {
          SetCellSizesAction(XCellSizes);
          XCellSizesChanged = true;
          return XCellSizes;
        }
      //! read-only access to the cell size in y-direction in m
      const t3DModelDim &GetYCellSizes() const
        {
          return YCellSizes;
        }
      //! read-write access to the cell size in y-direction in m
      t3DModelDim &SetYCellSizes()
        {
          SetCellSizesAction(YCellSizes);
          YCellSizesChanged = true;
          return YCellSizes;
        }
      //! read-only access to the cell size in z-direction in m
      const t3DModelDim &GetZCellSizes() const
        {
          return ZCellSizes;
        }
      //! read-write access to the cell size in z-direction in m
      t3DModelDim &SetZCellSizes()
        {
          SetCellSizesAction(ZCellSizes);
          ZCellSizesChanged = true;
          return ZCellSizes;
        }
      ThreeDModelBase();
      virtual ~ThreeDModelBase();
      };

  }

#endif /*THREEDMODELBASE_H_*/
