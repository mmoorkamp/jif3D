//============================================================================
// Name        : VTKTools.cpp
// Author      : Max Moorkamp
// Version     :
// Copyright   : 2008, MM
//============================================================================

#include "VTKTools.h"
#include <fstream>
#include "../Global/FatalException.h"

namespace jiba
  {

    //! A helper function to write the Coordinates for the 3 axes into a .vtk file
    void WriteCoordinatesToVTK(std::ofstream &file, const std::string &name,
        const ThreeDModelBase::t3DModelDim &CellSizes)
      {
        //write the coordinate name, the number of cell boundary values and its type
        file << name << " " << CellSizes.size() + 1 << " double" << std::endl;
        //our setup implies a (0,0,0) coordinate origin
        file << 0.0 << " ";
        //calculate the coordinates from the cell sizes and write to file
        std::partial_sum(CellSizes.begin(), CellSizes.end(),
            std::ostream_iterator<double>(file, " "));
        file << "\n";
      }

    void Write3DModelToVTK(const std::string &filename,
        const std::string &DataName,
        const ThreeDModelBase::t3DModelDim &XCellSizes,
        const ThreeDModelBase::t3DModelDim &YCellSizes,
        const ThreeDModelBase::t3DModelDim &ZCellSizes,
        const ThreeDModelBase::t3DModelData &Data)
      {
        //do some consistency checkes
        assert(Data.num_dimensions() == 3);
        assert(Data.shape()[0] == XCellSizes.size());
        assert(Data.shape()[1] == YCellSizes.size());
        assert(Data.shape()[2] == ZCellSizes.size());

        //get the size of each model direction
        const size_t nxvalues = XCellSizes.size();
        const size_t nyvalues = YCellSizes.size();
        const size_t nzvalues = ZCellSizes.size();

        std::ofstream outfile(filename.c_str());
        //first we have to write some general information about the file format
        outfile << "# vtk DataFile Version 2.0" << std::endl;
        outfile << "3D Model data" << std::endl;
        outfile << "ASCII" << std::endl;
        outfile << "DATASET RECTILINEAR_GRID" << std::endl;
        //we write the left and right boundaries of each cell, but only store the left
        //so we have to write one extra value
        outfile << "DIMENSIONS " << nxvalues + 1 << " " << nyvalues + 1 << " "
            << nzvalues + 1 << std::endl;
        //write information about the coordinate axes
        WriteCoordinatesToVTK(outfile, "X_COORDINATES", XCellSizes);
        WriteCoordinatesToVTK(outfile, "Y_COORDINATES", YCellSizes);
        WriteCoordinatesToVTK(outfile, "Z_COORDINATES", ZCellSizes);
        //write some information about the data itself
        outfile << "CELL_DATA " << nxvalues * nyvalues * nzvalues << std::endl;
        outfile << "SCALARS " << DataName << " double" << std::endl;
        outfile << "LOOKUP_TABLE default" << std::endl;
        //and then just the data values
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

    void Write3DDataToVTK(const std::string &filename,
        const std::string &DataName,
        const ThreeDGravityModel::tScalarMeasVec &Data,
        const ThreeDGravityModel::tMeasPosVec &PosX,
        const ThreeDGravityModel::tMeasPosVec &PosY,
        const ThreeDGravityModel::tMeasPosVec &PosZ)
      {
        //do some consistency checks
        const size_t ndata = Data.size();
        assert(ndata == PosX.size());
        assert(ndata == PosY.size());
        assert(ndata == PosZ.size());

        std::ofstream outfile(filename.c_str());
        //first we have to write some general information about the file format
        outfile << "# vtk DataFile Version 2.0" << std::endl;
        outfile << "3D Data" << std::endl;
        outfile << "ASCII" << std::endl;
        outfile << "DATASET UNSTRUCTURED_GRID" << std::endl;
        //we write the number of measurement points
        outfile << "POINTS " << ndata << " double" << std::endl;
        //write the coordinates of each point
        for (size_t i = 0; i < ndata; ++i)
          {
            outfile << PosX.at(i) << " " << PosY.at(i) << " " << PosZ.at(i)
                << "\n";
          }
        outfile << "POINT_DATA " << ndata << std::endl;
        outfile << "SCALARS " << DataName << " double" << std::endl;
        outfile << "LOOKUP_TABLE default" << std::endl;
        //and then just the data values
        std::copy(Data.begin(),Data.end(),std::ostream_iterator<double>(outfile," "));
        outfile << std::endl;
        if (outfile.fail())
          throw FatalException("Problem writing vtk  file");
      }

  }
