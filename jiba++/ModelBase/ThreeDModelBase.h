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
    /** \addtogroup modelbase Basic classes and routines for 3D models */
    /* @{ */

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
      //! The x-coordinate of the upper left front corner
      mutable t3DModelDim GridYCoordinates;
      //! The x-coordinate of the upper left front corner
      mutable t3DModelDim GridZCoordinates;
      //! Calculate the coordinates of the model cells from the sizes of each cell
      /*! This functions assumes that the coordinate of the upper left front corner of the model is
       * is (0,0,0).
       * @param Coordinates This vector will contain the coordinates of the left upper front corner of each cell
       * @param Sizes The size of each cell in m
       * @param ChangeFlag The flag that stores whether this coordinate has been changed
       */
      void CalcCoordinates(t3DModelDim &Coordinates, const t3DModelDim Sizes,
          bool &ChangeFlag) const
        {
          const size_t nelements = Sizes.size();
          if (ChangeFlag && nelements > 0)
            {
              Coordinates.resize(boost::extents[nelements]);
              if (nelements > 0)
                {
                  std::partial_sum(Sizes.begin(), Sizes.end(),
                      Coordinates.begin());
                  std::rotate(Coordinates.begin(), Coordinates.end() - 1,
                      Coordinates.end());
                  Coordinates[0] = 0.0;
                }
              ChangeFlag = false;
            }
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
      //!Get the x (north) coordinates of the cells, might perform calculations and write operation so it is not thread-safe
      const t3DModelDim &GetXCoordinates() const
        {
          CalcCoordinates(GridXCoordinates, XCellSizes, XCellSizesChanged);
          return GridXCoordinates;
        }
      //!Get the y (east) coordinates of the cells, might perform calculations and write operation so it is not thread-safe
      const t3DModelDim &GetYCoordinates() const
        {
          CalcCoordinates(GridYCoordinates, YCellSizes, YCellSizesChanged);
          return GridYCoordinates;
        }
      //!Get the z (depth) coordinates of the cells, might perform calculations and write operation so it is not thread-safe
      const t3DModelDim &GetZCoordinates() const
        {
          CalcCoordinates(GridZCoordinates, ZCellSizes, ZCellSizesChanged);
          return GridZCoordinates;
        }
      boost::array<ThreeDModelBase::t3DModelData::index,3>
          FindAssociatedIndices(const double xcoord, const double ycoord,
              const double zcoord) const;
      ThreeDModelBase();
      virtual ~ThreeDModelBase();
      };
  /* @} */
  }

#endif /*THREEDMODELBASE_H_*/
