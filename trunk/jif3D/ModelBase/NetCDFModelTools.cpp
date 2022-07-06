//============================================================================
// Name        : NetCDFModelTools.cpp
// Author      : Max Moorkamp
// Version     :
// Copyright   : 2008, MM
//============================================================================

#include "../Global/FatalException.h"
#include "../Global/NetCDFTools.h"
#include "../Global/NetCDFPortHelper.h"
#include "NetCDFModelTools.h"

#include <cassert>
#include <vector>
#include <algorithm>
using netCDF::NcDim;
using netCDF::NcFile;
using netCDF::NcVar;
using netCDF::NcVarAtt;

static const std::string OriginSuffix = "_Origin";

namespace jif3D
  {
    //This internal function should not appear in the doxygen documentation
    // Read a single CellCoordinates variable from a netcdf file
    /* This is a helper function that reads the cell coordinates for a single
     * dimension from the file.
     * @param NetCDFFile A netcdf file object ready for reading
     * @param CoordName The name of the dimension that specifies the cell sizes, e.g. "Northing"
     * @param CellCoord This will contain the cell sizes after the call
     * @return The number of cells
     */
    size_t ReadCoordinatesFromNetCDF(const NcFile &NetCDFFile,
        const std::string &CoordName, ThreeDModelBase::t3DModelDim &CellCoord)
      {
        //we store the information about the grid in a different way to
        //which we use it, we assume the first cell to have the upper left front corner at (0,0,0)
        //and store the right back lower coordinate for each cell
        //internally however we work with the upper left front corner and use the size
        //of the last cell to complete the geometry of the model

        //create a netcdf dimension with the chosen name
        NcDim CoordDim = NetCDFFile.getDim(CoordName.c_str());
        //determine the size of that dimension
        const size_t nvalues = CoordDim.getSize();
        //it does not make much sense to have an 0 size coordinate in a file
        //so we throw an exception
        if (nvalues == 0)
          {
            throw jif3D::FatalException("Cell coordinate is empty: " + CoordName,
            __FILE__,
            __LINE__);
          }
        //allocate memory in the class variable
        CellCoord.resize(nvalues + 1);

        // create netcdf variable with the same name as the dimension
        NcVar CoordVar = NetCDFFile.getVar(CoordName.c_str());

        //read coordinate values from netcdf file
        CoordVar.getVar(std::vector<std::size_t>(
          { 0 }), std::vector<std::size_t>(
          { nvalues }), CellCoord.data() + 1);

        double Origin = 0.0;
        std::string OriginName = CoordName + OriginSuffix;
        if (CheckExists(NetCDFFile, OriginName))
          {
            NcVar OVar = NetCDFFile.getVar(OriginName);
            OVar.getVar(&Origin);
          }
        CellCoord.at(0) = Origin;
//        SizeVar.get(&CellCoordinates[0], nvalues);
        //check whether the cell coordinates are sorted
        //otherwise we will have a problem
        if (!std::is_sorted(CellCoord.begin(), CellCoord.end()))
          {
            std::copy(CellCoord.begin(), CellCoord.end(),
                std::ostream_iterator<double>(std::cerr, " "));
            std::cerr << "\n";
            throw jif3D::FatalException(
                "Cell coordinates in netcdf file are not in increasing order: "
                    + CoordName, __FILE__,
                __LINE__);
          }

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
    NcDim WriteCoordinatesToNetCDF(NcFile &NetCDFFile, const std::string &CoordName,
        const ThreeDModelBase::t3DModelDim &CellCoord, const NcDim &BoundaryDim)
      {
        const size_t nvalues = CellCoord.size() - 1;
        // Add a dimension and a variable with the same name to the netcdf file
        NcDim CoordDim = NetCDFFile.addDim(CoordName, nvalues);
        NcVar CoordVar = NetCDFFile.addVar(CoordName, netCDF::ncDouble, CoordDim);
        //All length is measured in meters
        CoordVar.putAtt("units", "m");
        //We also store the name
        CoordVar.putAtt("long_name", (CoordName + " coordinate"));
        //we have to store the boundary values of the cells for plotting
        const std::string boundary_name = CoordName + "_bnds";
        CoordVar.putAtt("bounds", boundary_name);

        CoordVar.putVar(std::vector<std::size_t>(
          { 0 }), std::vector<std::size_t>(
          { nvalues }), CellCoord.data() + 1);

        //create the boundary variable
        std::vector<NcDim> dims;
        dims.push_back(CoordDim);
        dims.push_back(BoundaryDim);

        NcVar BoundaryVar = NetCDFFile.addVar(boundary_name, netCDF::ncDouble, dims);
        BoundaryVar.putAtt("units", "m");

        //the boundary variable contains two values per cell, the lower and upper boundary
        std::vector<double> BoundaryValues(nvalues * 2, 0);

        for (size_t i = 0; i < CellCoord.size() - 1; ++i)
          {
            BoundaryValues.at(2 * i) = CellCoord[i];
            BoundaryValues.at(2 * i + 1) = CellCoord[i + 1];
          }
//        BoundaryVar->put(&BoundaryValues[0], nvalues, 2); // old
        cxxport::put_legacy_ncvar(BoundaryVar, BoundaryValues.data(), nvalues, 2);

        std::string OriginName = CoordName + OriginSuffix;
        NcVar OrVar = NetCDFFile.addVar(OriginName, netCDF::ncDouble);
        OrVar.putVar(&CellCoord[0]);

        // We return the NcDim object, because we need it to write the model data
        return CoordDim;
      }

    /*! Read a  rectangular mesh  3D Model from a netcdf file
     * @param NetCDFFile An open NcFile object
     * @param DataName The name of the model data to retrieve, this allows several different models in one file
     * @param UnitsName The units of the model data, has to match the information in the file
     * @param XCellCoords The coordinates of the cells in x-direction in m (set by this routine)
     * @param YCellCoords The coordinates of the cells in y-direction in m (set by this routine)
     * @param ZCellCoords The coordinates of the cells in z-direction in m (set by this routine)
     * @param Data The model data (set by this routine)
     * @param xorigin The x-coordinate (Northing) of the (0,0,0) position of the mesh in m
     * @param yorigin The y-coordinate (Northing) of the (0,0,0) position of the mesh in m
     * @param zorigin The x-coordinate (Northing) of the (0,0,0) position of the mesh in m
     */
    void Read3DModelFromNetCDF(const NcFile &NetCDFFile, const std::string &DataName,
        const std::string &UnitsName, ThreeDModelBase::t3DModelDim &XCellCoords,
        ThreeDModelBase::t3DModelDim &YCellCoords,
        ThreeDModelBase::t3DModelDim &ZCellCoords, ThreeDModelBase::t3DModelData &Data)
      {

        //create netcdf variable for data
        NcVar DataVar = NetCDFFile.getVar(DataName);
        if (!DataVar.isNull())
          {
            //Read the sizes of the blocks in x,y and z-direction from the file
            const size_t nxvalues = ReadCoordinatesFromNetCDF(NetCDFFile, "Northing",
                XCellCoords);
            const size_t nyvalues = ReadCoordinatesFromNetCDF(NetCDFFile, "Easting",
                YCellCoords);
            const size_t nzvalues = ReadCoordinatesFromNetCDF(NetCDFFile, "Depth",
                ZCellCoords);
            //allocate memory for the data
            Data.resize(boost::extents[nxvalues][nyvalues][nzvalues]);
            //make sure we have the right units
            NcVarAtt Unit_Att = DataVar.getAtt("units");
            std::string unitInFile;
            Unit_Att.getValues(unitInFile);

            // if the units in the file are different from what we expect|
            if (unitInFile.compare(UnitsName) != 0)
              {
                throw jif3D::FatalException("Units in file do not match expected units !",
                    __FILE__,
                    __LINE__);
              }

            //read netcdf data from file
            //we store the data in the order Depth, Easting, Northing
            //but internally we work with x,y,z storage ordering
            std::vector<double> databuffer(nxvalues * nyvalues * nzvalues);

            cxxport::get_legacy_ncvar(DataVar, databuffer.data(), nzvalues, nyvalues,
                nxvalues);

            for (size_t i = 0; i < nzvalues; ++i)
              for (size_t j = 0; j < nyvalues; ++j)
                for (size_t k = 0; k < nxvalues; ++k)
                  Data[k][j][i] =
                      databuffer[k + j * nxvalues + i * (nxvalues * nyvalues)];

          }
        else
          {
            throw jif3D::FatalException("Cannot read variable: " + DataName, __FILE__,
            __LINE__);
          }
      }
    /*! Write the information about a rectangular mesh 3D model to a netcdf file
     * @param NetCDFFile An open NcFile object that allows writing
     * @param DataName The name of the model data (e.g. density, resistivity)
     * @param UnitsName The units of the model data
     * @param XCellCoords The coordinates of the cells in x-direction in m
     * @param YCellCoords The coordinates of the cells in y-direction in m
     * @param ZCellCoords The coordinates of the cells in z-direction in m
     * @param Data The model data, the shape has to match the length of the cell size vectors
     */
    void Write3DModelToNetCDF(NcFile &NetCDFFile, const std::string &DataName,
        const std::string &UnitsName, const ThreeDModelBase::t3DModelDim &XCellCoords,
        const ThreeDModelBase::t3DModelDim &YCellCoords,
        const ThreeDModelBase::t3DModelDim &ZCellCoords,
        const ThreeDModelBase::t3DModelData &Data)
      {
        //Make sure our data has the right dimensions and we have size information for each cell
        assert(Data.num_dimensions() == 3);
        assert(Data.shape()[0] == XCellCoords.size() - 1);
        assert(Data.shape()[1] == YCellCoords.size() - 1);
        assert(Data.shape()[2] == ZCellCoords.size() - 1);

        //add information about boundary values
        NcDim BoundaryDim = NetCDFFile.addDim("nbound", 2);
        NcVar BoundaryVar = NetCDFFile.addVar("nbound", netCDF::ncInt, BoundaryDim);
        BoundaryVar.putAtt("long_name", "Boundary index variable");

        std::vector<int> BIndices;
        BIndices.push_back(1);
        BIndices.push_back(2);
        cxxport::put_legacy_ncvar(BoundaryVar, BIndices.data(), BIndices.size());
//        BoundaryVar.put(&BIndices[0], 2); // old

        // Write the size information in x,y, and z-direction
        NcDim XCoordDim = WriteCoordinatesToNetCDF(NetCDFFile, "Northing", XCellCoords,
            BoundaryDim);
        NcDim YCoordDim = WriteCoordinatesToNetCDF(NetCDFFile, "Easting", YCellCoords,
            BoundaryDim);
        NcDim ZCoordDim = WriteCoordinatesToNetCDF(NetCDFFile, "Depth", ZCellCoords,
            BoundaryDim);

        std::vector<NcDim> dimVec;
        dimVec.push_back(ZCoordDim);
        dimVec.push_back(YCoordDim);
        dimVec.push_back(XCoordDim);

        NcVar DataVar = NetCDFFile.addVar(DataName, netCDF::ncDouble, dimVec);
        DataVar.putAtt("units", UnitsName);
        DataVar.putAtt("long_name", DataName);
        //Write the model data itself
        //we store the data in the order Depth, Easting, Northing
        //but internally we work with x,y,z storage ordering
        //so we allocate a temporary buffer
        std::vector<double> databuffer(Data.num_elements());
        //change the storage order
        const size_t xsize = XCoordDim.getSize();
        const size_t ysize = YCoordDim.getSize();
        const size_t zsize = ZCoordDim.getSize();

        for (size_t i = 0; i < zsize; ++i)
          for (size_t j = 0; j < ysize; ++j)
            for (size_t k = 0; k < xsize; ++k)
              databuffer[k + j * xsize + i * (xsize * ysize)] = Data[k][j][i];
        //write the buffer to the file
//        DataVar.put(databuffer, zsize, ysize, xsize);
        cxxport::put_legacy_ncvar(DataVar, databuffer.data(), zsize, ysize, xsize);

        NetCDFFile.putAtt("Conventions", "CF-1.3");
      }

    void ReadMeasPosNetCDF(const std::string &filename, std::vector<double> &PosX,
        std::vector<double> &PosY, std::vector<double> &PosZ)
      {

        //open the file
        NcFile DataFile(filename, NcFile::read);
        //read the three coordinates for the measurements
        ReadVec(DataFile, MeasPosXName, PosX);
        ReadVec(DataFile, MeasPosYName, PosY);
        ReadVec(DataFile, MeasPosZName, PosZ);
        //and make sure everything is consistent
        assert(PosX.size() == PosY.size());
        assert(PosX.size() == PosZ.size());
      }
  }
