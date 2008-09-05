//============================================================================
// Name        : NetCDFTools.cpp
// Author      : Max Moorkamp
// Version     :
// Copyright   : 2008, MM
//============================================================================

#include "NetCDFTools.h"

namespace jiba
  {

    size_t ReadSizesFromNetCDF(const NcFile &NetCDFFile,
        const std::string &SizeName, ThreeDModelBase::t3DModelDim &CellSize)
      {
        //we store the information about the grid in a different way to
        //which we use it, we assume the first cell to have the upper left front corner at (0,0,0)
        //and store the right back lower coordinate for each cell
        //internally however we work with the upper left front corner and use the size
        //of the last cell to complete the geometry of the model
        //Make sure we are reading into a one-dimensional entry
        assert(CellSize.num_dimensions() == 1);
        //create a netcdf dimension with the chosen name
        NcDim *SizeDim = NetCDFFile.get_dim(SizeName.c_str());
        //determine the size of that dimension
        const size_t nvalues = SizeDim->size();
        //allocate memory in the class variable
        CellSize.resize(boost::extents[nvalues]);
        ThreeDModelBase::t3DModelDim CellCoordinates(boost::extents[nvalues]);
        // create netcdf variable with the same name as the dimension
        NcVar *SizeVar = NetCDFFile.get_var(SizeName.c_str());
        //read coordinate values from netcdf file
        SizeVar->get(CellCoordinates.origin(), nvalues);
        // transform coordinates to cell sizes, assuming the start is at 0
        std::adjacent_difference(CellCoordinates.begin(),
            CellCoordinates.end(), CellSize.begin());
        //return the size of the current dimension
        return nvalues;
      }

    NcDim *WriteSizesToNetCDF(NcFile &NetCDFFile, const std::string &SizeName,
        const ThreeDModelBase::t3DModelDim &CellSize)
      {
        //Make sure we are writing a one-dimensional entry
        assert(CellSize.num_dimensions() == 1);
        // Add a dimension and a variable with the same name to the netcdf file
        NcDim *SizeDim = NetCDFFile.add_dim(SizeName.c_str(), CellSize.size());
        NcVar *SizeVar =
            NetCDFFile.add_var(SizeName.c_str(), ncDouble, SizeDim);
        //All length is measured in meters
        SizeVar->add_att("units", "m");
        //We also store the name
        SizeVar->add_att("long_name", (SizeName + " coordinate").c_str());
        // we store the coordinates of the cells in the netcdf file
        ThreeDModelBase::t3DModelDim CellCoordinates(
            boost::extents[CellSize.size()]);
        std::partial_sum(CellSize.begin(), CellSize.end(),
            CellCoordinates.begin());
        //Write the values
        SizeVar->put(CellCoordinates.origin(), CellCoordinates.size());
        // We return the NcDim object, because we need it to write the model data
        return SizeDim;
      }
    /*! Read a  rectangular mesh  3D Model from a netcdf file
     * @param NetCDFFile An open NcFile object
     * @param DataName The name of the model data to retrieve, this allows several different models in one file
     * @param UnitsName The units of the model data, has to match the information in the file
     * @param XCellSizes The sizes of the cells in x-direction in m (set by this routine)
     * @param YCellSizes The sizes of the cells in y-direction in m (set by this routine)
     * @param ZCellSizes The sizes of the cells in z-direction in m (set by this routine)
     * @param Data The model data (set by this routine)
     */
    void Read3DModelFromNetCDF(const NcFile &NetCDFFile,
        const std::string &DataName, const std::string &UnitsName,
        ThreeDModelBase::t3DModelDim &XCellSizes,
        ThreeDModelBase::t3DModelDim &YCellSizes,
        ThreeDModelBase::t3DModelDim &ZCellSizes,
        ThreeDModelBase::t3DModelData &Data)
      {

        //Read the sizes of the blocks in x,y and z-direction from the file
        const size_t nxvalues =
            ReadSizesFromNetCDF(NetCDFFile, "x", XCellSizes);
        const size_t nyvalues =
            ReadSizesFromNetCDF(NetCDFFile, "y", YCellSizes);
        const size_t nzvalues =
            ReadSizesFromNetCDF(NetCDFFile, "z", ZCellSizes);
        //allocate memory for the data
        Data.resize(boost::extents[nxvalues][nyvalues][nzvalues]);
        //create netcdf variable for data
        NcVar *DataVar = NetCDFFile.get_var(DataName.c_str());
        //make sure we have the right units
        NcAtt *Unit_Att = DataVar->get_att("units");
        std::string UnitInFile = Unit_Att->as_string(0);
        // if the units in the file are different from what we expect|
        if (UnitInFile.compare(UnitsName) != 0)
          throw std::runtime_error(
              "Units in file do not match expected units !");
        //read netcdf data from file
        DataVar->get(Data.origin(), nxvalues, nyvalues, nzvalues);
      }
    /*! Write the information about a rectangular mesh 3D model to a netcdf file
     * @param NetCDFFile An open NcFile object that allows writing
     * @param DataName The name of the model data (e.g. density, resistivity)
     * @param UnitsName The units of the model data
     * @param XCellSizes The sizes of the cells in x-direction in m
     * @param YCellSizes The sizes of the cells in y-direction in m
     * @param ZCellSizes The sizes of the cells in z-direction in m
     * @param Data The model data, the shape has to match the length of the cell size vectors
     */
    void Write3DModelToNetCDF(NcFile &NetCDFFile, const std::string &DataName,
        const std::string &UnitsName,const ThreeDModelBase::t3DModelDim &XCellSizes,
        const ThreeDModelBase::t3DModelDim &YCellSizes,
        const ThreeDModelBase::t3DModelDim &ZCellSizes,
        const ThreeDModelBase::t3DModelData &Data)
      {
        //Make sure our data has the right dimensions and we have size information for each cell
        assert(Data.num_dimensions() == 3);
        assert(Data.shape()[0] == XCellSizes.size());
        assert(Data.shape()[1] == YCellSizes.size());
        assert(Data.shape()[2] == ZCellSizes.size());
        // Write the size information in x,y, and z-direction
        NcDim *XSizeDim = WriteSizesToNetCDF(NetCDFFile, "x", XCellSizes);
        NcDim *YSizeDim = WriteSizesToNetCDF(NetCDFFile, "y", YCellSizes);
        NcDim *ZSizeDim = WriteSizesToNetCDF(NetCDFFile, "z", ZCellSizes);
        NcVar *DataVar = NetCDFFile.add_var(DataName.c_str(), ncDouble,
            XSizeDim, YSizeDim, ZSizeDim);
        DataVar->add_att("units", UnitsName.c_str());
        //Write the model data itself
        DataVar->put(Data.origin(), XSizeDim->size(), YSizeDim->size(),
            ZSizeDim->size());
      }
  }
