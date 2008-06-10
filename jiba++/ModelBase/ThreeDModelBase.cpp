//============================================================================
// Name        : ThreeDModelBase.cpp
// Author      : Max Moorkamp
// Version     :
// Copyright   : 2008, MM
//============================================================================

#include "ThreeDModelBase.h"
#include <cassert>
#include <fstream>
#include "../Global/FatalException.h"

namespace jiba
  {

    ThreeDModelBase::ThreeDModelBase() :
      XCellSizesChanged(true), YCellSizesChanged(true), ZCellSizesChanged(true)
      {
      }

    ThreeDModelBase::~ThreeDModelBase()
      {
      }

    size_t ThreeDModelBase::ReadSizesFromNetCDF(const NcFile &NetCDFFile,
        const std::string &SizeName, t3DModelDim &CellSize,
        t3DModelDim CellCoordinates)
      {
        //Make sure we are reading into a one-dimensional entry
        assert(CellSize.num_dimensions() == 1);
        //create a netcdf dimension with the chosen name
        NcDim *SizeDim = NetCDFFile.get_dim(SizeName.c_str());
        //determine the size of that dimension
        const size_t nvalues = SizeDim->size();
        //allocate memory in the class variable
        CellSize.resize(boost::extents[nvalues]);
        CellCoordinates.resize(boost::extents[nvalues]);
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

    NcDim *ThreeDModelBase::WriteSizesToNetCDF(NcFile &NetCDFFile,
        const std::string &SizeName, const t3DModelDim &CellSize) const
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
        SizeVar->add_att("long_name", (SizeName+" coordinate").c_str());
        // we store the coordinates of the cells in the netcdf file
        t3DModelDim CellCoordinates(boost::extents[CellSize.size()]);
        std::partial_sum(CellSize.begin(), CellSize.end(),
            CellCoordinates.begin());
        //Write the values
        SizeVar->put(CellCoordinates.origin(), CellCoordinates.size());
        // We return the NcDim object, because we need it to write the model data
        return SizeDim;
      }

    void ThreeDModelBase::ReadDataFromNetCDF(const NcFile &NetCDFFile,
        const std::string &DataName, const std::string &UnitsName)
      {

        //Read the sizes of the blocks in x,y and z-direction from the file
        const size_t nxvalues = ReadSizesFromNetCDF(NetCDFFile, "x",
            XCellSizes, GridXCoordinates);
        const size_t nyvalues = ReadSizesFromNetCDF(NetCDFFile, "y",
            YCellSizes, GridYCoordinates);
        const size_t nzvalues = ReadSizesFromNetCDF(NetCDFFile, "z",
            ZCellSizes, GridZCoordinates);
        //allocate memory for the data
        Data.resize(boost::extents[nxvalues][nyvalues][nzvalues]);
        //create netcdf variable for data
        NcVar *DataVar = NetCDFFile.get_var(DataName.c_str());
        //make sure we have the right units
        NcAtt *Unit_Att = DataVar->get_att("units");
        std::string UnitInFile = Unit_Att->as_string(0);
        // if the units in the file are different from what we expect|
        if (UnitInFile.compare(UnitsName) != 0)
          throw std::runtime_error("Units in file do not match expected units !");
        //read netcdf data from file
        DataVar->get(Data.origin(), nxvalues, nyvalues, nzvalues);
      }

    void ThreeDModelBase::WriteDataToNetCDF(NcFile &NetCDFFile,
        const std::string &DataName, const std::string &UnitsName) const
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

    void ThreeDModelBase::WriteVTK(std::string filename)
      {
        assert(Data.num_dimensions() == 3);
        assert(Data.shape()[0] == XCellSizes.size());
        assert(Data.shape()[1] == YCellSizes.size());
        assert(Data.shape()[2] == ZCellSizes.size());

        const size_t nxvalues = XCellSizes.size();
        const size_t nyvalues = YCellSizes.size();
        const size_t nzvalues = ZCellSizes.size();

        std::ofstream outfile(filename.c_str());
        outfile << "# vtk DataFile Version 2.0" << std::endl;
        outfile << "3D Model data" << std::endl;
        outfile << "ASCII" << std::endl;
        outfile << "DATASET RECTILINEAR_GRID" << std::endl;
        outfile << "DIMENSIONS " << nxvalues+1 << " " <<nyvalues+1 << " "
            << nzvalues +1 << std::endl;

        outfile << "X_COORDINATES "<< nxvalues+1 << " double" << std::endl;
        outfile << 0.0 << " ";
        std::partial_sum(GetXCellSizes().begin(), GetXCellSizes().end(), std::ostream_iterator<double>(outfile, " "));
        outfile << std::endl;

        outfile << "Y_COORDINATES "<< nyvalues+1 << " double" << std::endl;
        outfile << 0.0 << " ";
        std::partial_sum(GetYCellSizes().begin(), GetYCellSizes().end(), std::ostream_iterator<double>(outfile, " "));
        outfile << std::endl;

        outfile << "Z_COORDINATES "<< nzvalues+1 << " double" << std::endl;
        outfile << 0.0 << " ";
        std::partial_sum(GetZCellSizes().begin(), GetZCellSizes().end(), std::ostream_iterator<double>(outfile, " "));
        outfile << std::endl;

        
        outfile << "CELL_DATA " << nxvalues * nyvalues * nzvalues
            << std::endl;
        outfile << "SCALARS Resistivity double" << std::endl;
        outfile << "LOOKUP_TABLE default" << std::endl;
        for (size_t i = 0; i < nzvalues; ++i)
          {
            for (size_t j = 0; j < nyvalues; ++j)
              {
                for (size_t k = 0; k < nxvalues; ++k)
                  {
                    outfile << Data[k][j][i] << " ";
                  }
                outfile << std::endl;
              }
          }
        if (outfile.fail())
          throw FatalException("Problem writing vtk  file");
      }

  }
