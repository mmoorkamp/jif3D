//============================================================================
// Name        : NetCDFTools.cpp
// Author      : Max Moorkamp
// Version     :
// Copyright   : 2008, MM
//============================================================================

#include "NetCDFTools.h"
#include <boost/assign/std/vector.hpp>

namespace jiba
  {
    //This internal function should not appear in the doxygen documentation
    // Read a single CellSize variable from a netcdf file
    /* This is a helper function that reads the cell length for a single
     * dimension from the file.
     * @param NetCDFFile A netcdf file object ready for reading
     * @param SizeName The name of the dimension that specifies the cell sizes, e.g. "Northing"
     * @param CellSize This will contain the cell sizes after the call
     * @return The number of cells
     */
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
    //This internal function should not appear in the doxygen documentation
    // Write the length of the model cell along a single dimension to the file
    /* A helper function to write the cell sizes and the corresponding
     * dimension to a netcdf file.
     * @param NetCDFFile A netcdf file object ready for writing
     * @param SizeName  The name of the cell size variable, e.g. "Northing"
     * @param CellSize The vector of cell sizes
     * @param BoundaryDim For plotting we also create a boundary variable for each cell size variable, as the dimension for this is shared by all cell sizes, we create it once and pass it to this function
     * @return A pointer to the netcdf dimension object for the cell size variable
     */
    NcDim *WriteSizesToNetCDF(NcFile &NetCDFFile, const std::string &SizeName,
        const ThreeDModelBase::t3DModelDim &CellSize, NcDim *BoundaryDim)
      {
        //Make sure we are writing a one-dimensional entry
        assert(CellSize.num_dimensions() == 1);
        const size_t nvalues = CellSize.size();
        // Add a dimension and a variable with the same name to the netcdf file
        NcDim *SizeDim = NetCDFFile.add_dim(SizeName.c_str(), nvalues);
        NcVar *SizeVar =
            NetCDFFile.add_var(SizeName.c_str(), ncDouble, SizeDim);
        //All length is measured in meters
        SizeVar->add_att("units", "m");
        //We also store the name
        SizeVar->add_att("long_name", (SizeName + " coordinate").c_str());
        //we have to store the boundary values of the cells for plotting
        const std::string boundary_name = SizeName + "_bnds";
        SizeVar->add_att("bounds", boundary_name.c_str());
        // we store the coordinates of the cells in the netcdf file
        ThreeDModelBase::t3DModelDim CellCoordinates(boost::extents[nvalues]);
        std::partial_sum(CellSize.begin(), CellSize.end(),
            CellCoordinates.begin());
        //Write the values
        SizeVar->put(CellCoordinates.origin(), nvalues);
        //create the boundary variable
        NcVar *BoundaryVar = NetCDFFile.add_var(boundary_name.c_str(),
            ncDouble, SizeDim, BoundaryDim);
        BoundaryVar->add_att("units", "m");

        //the boundary variable contains two values per cell, the lower and upper boundary
        std::vector<double> BoundaryValues(nvalues * 2, 0);
        BoundaryValues.at(0) = 0.0;
        BoundaryValues.at(1) = CellCoordinates[0];

        for (size_t i = 1; i < CellCoordinates.size(); ++i)
          {
            BoundaryValues.at(2 * i) = CellCoordinates[i - 1];
            BoundaryValues.at(2 * i + 1) = CellCoordinates[i];
          }
        BoundaryVar->put(BoundaryValues.data(), nvalues, 2);
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
        const size_t nxvalues = ReadSizesFromNetCDF(NetCDFFile, "Northing",
            XCellSizes);
        const size_t nyvalues = ReadSizesFromNetCDF(NetCDFFile, "Easting",
            YCellSizes);
        const size_t nzvalues = ReadSizesFromNetCDF(NetCDFFile, "Depth",
            ZCellSizes);
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
        //we store the data in the order Depth, Easting, Northing
        //but internally we work with x,y,z storage ordering
        double *databuffer = new double[nxvalues * nyvalues * nzvalues];

        DataVar->get(databuffer, nzvalues, nyvalues, nxvalues);
        for (size_t i = 0; i < nzvalues; ++i)
          for (size_t j = 0; j < nyvalues; ++j)
            for (size_t k = 0; k < nxvalues; ++k)
              Data[k][j][i] = databuffer[k + j * nxvalues + i * (nxvalues
                  * nyvalues)];
        delete[] databuffer;
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
        const std::string &UnitsName,
        const ThreeDModelBase::t3DModelDim &XCellSizes,
        const ThreeDModelBase::t3DModelDim &YCellSizes,
        const ThreeDModelBase::t3DModelDim &ZCellSizes,
        const ThreeDModelBase::t3DModelData &Data)
      {
        using namespace boost::assign;
        //Make sure our data has the right dimensions and we have size information for each cell
        assert(Data.num_dimensions() == 3);
        assert(Data.shape()[0] == XCellSizes.size());
        assert(Data.shape()[1] == YCellSizes.size());
        assert(Data.shape()[2] == ZCellSizes.size());

        //add information about boundary values
        NcDim *BoundaryDim = NetCDFFile.add_dim("nbound", 2);
        NcVar *BoundaryVar = NetCDFFile.add_var("nbound", ncInt, BoundaryDim);
        BoundaryVar->add_att("long_name", "Boundary index variable");
        std::vector<int> BIndices;
        BIndices += 1, 2;
        BoundaryVar->put(&BIndices[0], 2);
        // Write the size information in x,y, and z-direction
        NcDim *XSizeDim = WriteSizesToNetCDF(NetCDFFile, "Northing",
            XCellSizes, BoundaryDim);
        NcDim *YSizeDim = WriteSizesToNetCDF(NetCDFFile, "Easting", YCellSizes,
            BoundaryDim);
        NcDim *ZSizeDim = WriteSizesToNetCDF(NetCDFFile, "Depth", ZCellSizes,
            BoundaryDim);

        NcVar *DataVar = NetCDFFile.add_var(DataName.c_str(), ncDouble,
            ZSizeDim, YSizeDim, XSizeDim);
        DataVar->add_att("units", UnitsName.c_str());
        DataVar->add_att("long_name", DataName.c_str());
        //Write the model data itself
        //we store the data in the order Depth, Easting, Northing
        //but internally we work with x,y,z storage ordering
        //so we allocate a temporary buffer
        double *databuffer = new double[Data.num_elements()];
        //change the storage order
        const size_t xsize = XSizeDim->size();
        const size_t ysize = YSizeDim->size();
        const size_t zsize = ZSizeDim->size();
        for (size_t i = 0; i < zsize; ++i)
          for (size_t j = 0; j < ysize; ++j)
            for (size_t k = 0; k < xsize; ++k)
              databuffer[k + j * xsize + i * (xsize * ysize)] = Data[k][j][i];
        //write the buffer to the file
        DataVar->put(databuffer, ZSizeDim->size(), YSizeDim->size(),
            XSizeDim->size());
        //and delete it
        delete[] databuffer;
        NetCDFFile.add_att("Conventions", "CF-1.3");
      }

    void ReadMeasPosNetCDF(const std::string filename,
        ThreeDModelBase::tMeasPosVec &PosX, ThreeDModelBase::tMeasPosVec &PosY,
        ThreeDModelBase::tMeasPosVec &PosZ)
      {

        //open the file
        NcFile DataFile(filename.c_str(), NcFile::ReadOnly);
        //read the three coordinates for the measurements
        ReadVec(DataFile, MeasPosXName, StationNumberName,PosX);
        ReadVec(DataFile, MeasPosYName, StationNumberName,PosY);
        ReadVec(DataFile, MeasPosZName, StationNumberName,PosZ);
        //and make sure everything is consistent
        assert(PosX.size() == PosY.size());
        assert(PosX.size() == PosZ.size());
      }
  }
