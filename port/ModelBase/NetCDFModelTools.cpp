//============================================================================
// Name        : NetCDFModelTools.cpp
// Author      : Max Moorkamp
// Version     :
// Copyright   : 2008, MM
//============================================================================

#include <cassert>

#include "../Global/FatalException.h"
#include "../Global/NetCDFTools.h"
#include "../Global/NetCDFPortHelper.h"

#include "NetCDFModelTools.h"

using netCDF::NcDim;
using netCDF::NcFile;
using netCDF::NcVar;
using netCDF::NcVarAtt;

namespace jif3D
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
    size_t ReadSizesFromNetCDF(const NcFile &NetCDFFile, const std::string &SizeName,
        ThreeDModelBase::t3DModelDim &CellSize)
      {
        //we store the information about the grid in a different way to
        //which we use it, we assume the first cell to have the upper left front corner at (0,0,0)
        //and store the right back lower coordinate for each cell
        //internally however we work with the upper left front corner and use the size
        //of the last cell to complete the geometry of the model

        //create a netcdf dimension with the chosen name
        NcDim SizeDim = NetCDFFile.getDim(SizeName.c_str());
        //determine the size of that dimension
        const size_t nvalues = SizeDim.getSize();
        //it does not make much sense to have an 0 size coordinate in a file
        //so we throw an exception
        if (nvalues == 0)
          {
            throw jif3D::FatalException("Cell coordinate is empty: " + SizeName, __FILE__,
                __LINE__);
          }
        //allocate memory in the class variable
        CellSize.resize(nvalues);
        ThreeDModelBase::t3DModelDim CellCoordinates(nvalues);
        // create netcdf variable with the same name as the dimension
        NcVar SizeVar = NetCDFFile.getVar(SizeName.c_str());

        std::vector<size_t> starts(nvalues, 0u);

        //read coordinate values from netcdf file
        jif3D::cxxport::get_legacy_ncvar(SizeVar, CellCoordinates.data(), nvalues);
//        SizeVar.get(&CellCoordinates[0], nvalues);
        //check whether the cell coordinates are sorted
        //otherwise we will have a problem
        ThreeDModelBase::t3DModelDim::iterator pos = std::adjacent_find(
            CellCoordinates.begin(), CellCoordinates.end(), // range
            std::greater<int>());

        if (pos != CellCoordinates.end())
          {
            throw jif3D::FatalException(
                "Cell coordinates in netcdf file are not in increasing order.", __FILE__,
                __LINE__);
          }

        // transform coordinates to cell sizes, assuming the start is at 0
        std::adjacent_difference(CellCoordinates.begin(), CellCoordinates.end(),
            CellSize.begin());
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
    NcDim WriteSizesToNetCDF(NcFile &NetCDFFile, const std::string &SizeName,
        const ThreeDModelBase::t3DModelDim &CellSize, const NcDim &BoundaryDim)
      {
        const size_t nvalues = CellSize.size();
        // Add a dimension and a variable with the same name to the netcdf file
        NcDim SizeDim = NetCDFFile.addDim(SizeName, nvalues);
        NcVar SizeVar = NetCDFFile.addVar(SizeName, netCDF::ncDouble, SizeDim);
        //All length is measured in meters
        SizeVar.putAtt("units", "m");
        //We also store the name
        SizeVar.putAtt("long_name", (SizeName + " coordinate"));
        //we have to store the boundary values of the cells for plotting
        const std::string boundary_name = SizeName + "_bnds";
        SizeVar.putAtt("bounds", boundary_name);

        // we store the coordinates of the cells in the netcdf file
        ThreeDModelBase::t3DModelDim CellCoordinates(nvalues);
        std::partial_sum(CellSize.begin(), CellSize.end(), CellCoordinates.begin());
        //Write the values
//        SizeVar->put(&CellCoordinates[0], nvalues); // old way
        cxxport::put_legacy_ncvar(SizeVar, CellCoordinates.data(), nvalues);

        //create the boundary variable
        std::vector<NcDim> dims;
        dims.push_back(SizeDim);
        dims.push_back(BoundaryDim);

        NcVar BoundaryVar = NetCDFFile.addVar(boundary_name, netCDF::ncDouble, dims);
        BoundaryVar.putAtt("units", "m");

        //the boundary variable contains two values per cell, the lower and upper boundary
        std::vector<double> BoundaryValues(nvalues * 2, 0);
        BoundaryValues[0] = 0.0;
        BoundaryValues[1] = CellCoordinates[0];

        for (size_t i = 1u; i < CellCoordinates.size(); ++i)
          {
            BoundaryValues.at(2 * i) = CellCoordinates[i - 1];
            BoundaryValues.at(2 * i + 1) = CellCoordinates[i];
          }
//        BoundaryVar->put(&BoundaryValues[0], nvalues, 2); // old
        cxxport::put_legacy_ncvar(BoundaryVar, BoundaryValues.data(), nvalues, 2);
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
    void Read3DModelFromNetCDF(const NcFile &NetCDFFile, const std::string &DataName,
        const std::string &UnitsName, ThreeDModelBase::t3DModelDim &XCellSizes,
        ThreeDModelBase::t3DModelDim &YCellSizes,
        ThreeDModelBase::t3DModelDim &ZCellSizes, ThreeDModelBase::t3DModelData &Data)
      {

        //create netcdf variable for data
        NcVar DataVar = NetCDFFile.getVar(DataName);
        if (!DataVar.isNull())
          {
            //Read the sizes of the blocks in x,y and z-direction from the file
            const size_t nxvalues = ReadSizesFromNetCDF(NetCDFFile, "Northing",
                XCellSizes);
            const size_t nyvalues = ReadSizesFromNetCDF(NetCDFFile, "Easting",
                YCellSizes);
            const size_t nzvalues = ReadSizesFromNetCDF(NetCDFFile, "Depth", ZCellSizes);
            //allocate memory for the data
            Data.resize(boost::extents[nxvalues][nyvalues][nzvalues]);
            //make sure we have the right units
            NcVarAtt Unit_Att = DataVar.getAtt("units");
            std::string unitInFile;
            Unit_Att.getValues(unitInFile);

            // if the units in the file are different from what we expect|
            if (unitInFile.compare(UnitsName) != 0) {
              throw std::runtime_error("Units in file do not match expected units !");
            }

            //read netcdf data from file
            //we store the data in the order Depth, Easting, Northing
            //but internally we work with x,y,z storage ordering
            double *databuffer = new double[nxvalues * nyvalues * nzvalues];;

//            DataVar.get(databuffer, nzvalues, nyvalues, nxvalues);
            cxxport::get_legacy_ncvar(DataVar, databuffer, nzvalues, nyvalues, nxvalues);

            for (size_t i = 0; i < nzvalues; ++i)
              for (size_t j = 0; j < nyvalues; ++j)
                for (size_t k = 0; k < nxvalues; ++k)
                  Data[k][j][i] =
                      databuffer[k + j * nxvalues + i * (nxvalues * nyvalues)];

            delete[] databuffer;
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
     * @param XCellSizes The sizes of the cells in x-direction in m
     * @param YCellSizes The sizes of the cells in y-direction in m
     * @param ZCellSizes The sizes of the cells in z-direction in m
     * @param Data The model data, the shape has to match the length of the cell size vectors
     */
    void Write3DModelToNetCDF(NcFile &NetCDFFile, const std::string &DataName,
        const std::string &UnitsName, const ThreeDModelBase::t3DModelDim &XCellSizes,
        const ThreeDModelBase::t3DModelDim &YCellSizes,
        const ThreeDModelBase::t3DModelDim &ZCellSizes,
        const ThreeDModelBase::t3DModelData &Data)
      {
        //Make sure our data has the right dimensions and we have size information for each cell
        assert(Data.num_dimensions() == 3);
        assert(Data.shape()[0] == XCellSizes.size());
        assert(Data.shape()[1] == YCellSizes.size());
        assert(Data.shape()[2] == ZCellSizes.size());

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
        NcDim XSizeDim = WriteSizesToNetCDF(NetCDFFile, "Northing", XCellSizes,
            BoundaryDim);
        NcDim YSizeDim = WriteSizesToNetCDF(NetCDFFile, "Easting", YCellSizes,
            BoundaryDim);
        NcDim ZSizeDim = WriteSizesToNetCDF(NetCDFFile, "Depth", ZCellSizes,
            BoundaryDim);

        std::vector<NcDim> dimVec;
        dimVec.push_back(ZSizeDim);
        dimVec.push_back(YSizeDim);
        dimVec.push_back(XSizeDim);

        NcVar DataVar = NetCDFFile.addVar(DataName, netCDF::ncDouble, dimVec);
        DataVar.putAtt("units", UnitsName);
        DataVar.putAtt("long_name", DataName);
        //Write the model data itself
        //we store the data in the order Depth, Easting, Northing
        //but internally we work with x,y,z storage ordering
        //so we allocate a temporary buffer
        double *databuffer = new double[Data.num_elements()];
        //change the storage order
        const size_t xsize = XSizeDim.getSize();
        const size_t ysize = YSizeDim.getSize();
        const size_t zsize = ZSizeDim.getSize();

        for (size_t i = 0; i < zsize; ++i)
          for (size_t j = 0; j < ysize; ++j)
            for (size_t k = 0; k < xsize; ++k)
              databuffer[k + j * xsize + i * (xsize * ysize)] = Data[k][j][i];
        //write the buffer to the file
//        DataVar.put(databuffer, zsize, ysize, xsize);
        cxxport::put_legacy_ncvar(DataVar, databuffer, zsize, ysize, xsize);

        //and delete it
        delete[] databuffer;
        NetCDFFile.putAtt("Conventions", "CF-1.3");
      }

    void ReadMeasPosNetCDF(const std::string &filename, ThreeDModelBase::tMeasPosVec &PosX,
        ThreeDModelBase::tMeasPosVec &PosY, ThreeDModelBase::tMeasPosVec &PosZ)
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
