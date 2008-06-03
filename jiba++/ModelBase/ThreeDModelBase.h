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
      virtual void SetCellSizesAction(t3DModelDim &sizes){}
      bool XCellSizesChanged;
      bool YCellSizesChanged;
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
  protected:
      //The access to the data is protected, this requires the derived class
      // to provide an access function and allows to add initial checks and 
      // meaningfull names
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
      const t3DModelDim &GetXCoordinates()
        {
          if (XCellSizesChanged && XCellSizes.size() > 0)
            {
              GridXCoordinates.resize(boost::extents[XCellSizes.size()]);
              std::partial_sum(XCellSizes.begin(), XCellSizes.end(),
                  GridXCoordinates.begin());
              std::rotate(GridXCoordinates.begin(), GridXCoordinates.end()-1,
                  GridXCoordinates.end());
              GridXCoordinates[0] = 0.0;
              XCellSizesChanged = false;
            }
          return GridXCoordinates;
        }
      const t3DModelDim &GetYCoordinates()
        {
          if (YCellSizesChanged && YCellSizes.size() > 0)
            {
              GridYCoordinates.resize(boost::extents[YCellSizes.size()]);
              std::partial_sum(YCellSizes.begin(), YCellSizes.end(),
                  GridYCoordinates.begin());
              std::rotate(GridYCoordinates.begin(), GridYCoordinates.end()-1,
                  GridYCoordinates.end());
              GridYCoordinates[0] = 0.0;
              YCellSizesChanged = false;
            }
          return GridYCoordinates;
        }
      const t3DModelDim &GetZCoordinates()
        {
          if (ZCellSizesChanged && ZCellSizes.size() > 0)
            {
              GridZCoordinates.resize(boost::extents[ZCellSizes.size()]);
              std::partial_sum(ZCellSizes.begin(), ZCellSizes.end(),
                  GridZCoordinates.begin());
              std::rotate(GridZCoordinates.begin(), GridZCoordinates.end()-1,
                  GridZCoordinates.end());
              GridZCoordinates[0] = 0.0;
            }
          return GridZCoordinates;
        }
      //! Read data and associated cell sizes from a netcdf file
      void ReadDataFromNetCDF(const NcFile &NetCDFFile,
          const std::string &DataName, const std::string &UnitsName);
      //! Write data and associated cell sizes to a netcdf file
      void WriteDataToNetCDF(NcFile &NetCDFFile, const std::string &DataName,
          const std::string &UnitsName) const;
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
