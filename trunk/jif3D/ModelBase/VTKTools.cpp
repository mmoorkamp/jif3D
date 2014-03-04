//============================================================================
// Name        : VTKTools.cpp
// Author      : Max Moorkamp
// Version     :
// Copyright   : 2008, MM
//============================================================================

#include "VTKTools.h"
#include <fstream>
#include "../Global/FatalException.h"

 #pragma GCC diagnostic ignored "-Wuninitialized"

namespace jif3D
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

    //! A helper function to read the Coordinates for the 3 axes from a .vtk file
    void ReadCoordinatesTFromVTK(std::ifstream &file,
        ThreeDModelBase::t3DModelDim &CellSizes)
      {
        std::string dummy;
        int nvalues;
        file >> dummy >> nvalues >> dummy;
        std::vector<double> coordinates(nvalues);
        for (int i = 0; i < nvalues; ++i)
          {
            file >> coordinates[i];
          }
        CellSizes.resize(boost::extents[nvalues - 1]);
        std::adjacent_difference(coordinates.begin() + 1, coordinates.end(),
            CellSizes.begin());
      }

    /*! Write a .vtk file to plot a 3D model
     * @param filename The name of the output file, should contain the ending .vtk
     * @param DataName The name of the model data, for information for plotting programs
     * @param XCellSizes The sizes of the cells in x-direction in m
     * @param YCellSizes The sizes of the cells in y-direction in m
     * @param ZCellSizes The sizes of the cells in z-direction in m
     * @param Data The model values within each cell, shape has to match the  cell sizes
     */
    void Write3DModelToVTK(const std::string &filename, const std::string &DataName,
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

        typedef boost::multi_array_types::index_range range;
        ThreeDModelBase::t3DModelData::index_gen indices;
        for (size_t i = 0; i < nzvalues; ++i)
          {
            for (size_t j = 0; j < nyvalues; ++j)
              {
                ThreeDModelBase::t3DModelData::const_array_view<1>::type myview =
                    Data[indices[range(0, nxvalues)][j][i]];
                std::copy(myview.begin(), myview.end(),
                    std::ostream_iterator<double>(outfile, " "));
                outfile << "\n";
              }
          }
        if (outfile.fail())
          throw FatalException("Problem writing vtk  file: " + filename);
      }

    /*! Read a .vtk file containing a 3D model
     * @param filename The name of the input file, should contain the ending .vtk
     * @param XCellSizes The sizes of the cells in x-direction in m
     * @param YCellSizes The sizes of the cells in y-direction in m
     * @param ZCellSizes The sizes of the cells in z-direction in m
     * @param Data The model values within each cell, shape has to match the  cell sizes
     */
    void Read3DModelFromVTK(const std::string &filename,
        ThreeDModelBase::t3DModelDim &XCellSizes,
        ThreeDModelBase::t3DModelDim &YCellSizes,
        ThreeDModelBase::t3DModelDim &ZCellSizes, ThreeDModelBase::t3DModelData &Data)
      {
        std::ifstream infile(filename.c_str());
        char dummy[255];
        //read the header information, but we do not need it
        infile.getline(dummy, 255);
        infile.getline(dummy, 255);
        infile.getline(dummy, 255);
        infile.getline(dummy, 255);
        //get the number of values in each direction
        int nxvalues, nyvalues, nzvalues;
        infile >> dummy >> nxvalues >> nyvalues >> nzvalues;
        //for vtk we need both boundaries of the cells,
        //but internally we just use the right boundary
        //so we need on element less for each axis
        nxvalues -= 1;
        nyvalues -= 1;
        nzvalues -= 1;
        //read the coordinate of the axes
        ReadCoordinatesTFromVTK(infile, XCellSizes);
        ReadCoordinatesTFromVTK(infile, YCellSizes);
        ReadCoordinatesTFromVTK(infile, ZCellSizes);
        //skip some more vtk specific information
        infile.getline(dummy, 255);
        infile.getline(dummy, 255);
        infile.getline(dummy, 255);
        infile.getline(dummy, 255);
        Data.resize(boost::extents[nxvalues][nyvalues][nzvalues]);
        //read the values from the file, the storage ordering
        //is different from our internal storage, so we
        //have to use nested loops
        for (int i = 0; i < nzvalues; ++i)
          {
            for (int j = 0; j < nyvalues; ++j)
              {
                for (int k = 0; k < nxvalues; ++k)
                  {
                    infile >> Data[k][j][i];
                    if (infile.fail())
                      throw FatalException("Problem reading vtk  file: " + filename);
                  }
              }
          }
      }

    //helper function that writes the common header information
    //for scalar and tensor data
    void WriteDataHeader(std::ofstream &outfile, const ThreeDModelBase::tMeasPosVec &PosX,
        const ThreeDModelBase::tMeasPosVec &PosY,
        const ThreeDModelBase::tMeasPosVec &PosZ)
      {
        const size_t ndata = PosX.size();
        outfile << "# vtk DataFile Version 2.0" << std::endl;
        outfile << "3D Data" << std::endl;
        outfile << "ASCII" << std::endl;
        outfile << "DATASET UNSTRUCTURED_GRID" << std::endl;
        //we write the number of measurement points
        outfile << "POINTS " << ndata << " double" << std::endl;
        //write the coordinates of each point
        for (size_t i = 0; i < ndata; ++i)
          {
            outfile << PosX.at(i) << " " << PosY.at(i) << " " << PosZ.at(i) << "\n";
          }
        outfile << "POINT_DATA " << ndata << std::endl;
      }

    /*! Write a collection of scalar measurements to a .vtk file for plotting
     * @param filename The name of the file, should contain the ending .vtk
     * @param DataName The name of the data for information in the plotting program
     * @param Data The vector of data points, 1 datum per position
     * @param PosX The position of the measurement points in x-direction in m
     * @param PosY The position of the measurement points in y-direction in m
     * @param PosZ The position of the measurement points in z-direction in m
     */
    void Write3DDataToVTK(const std::string &filename, const std::string &DataName,
        const jif3D::rvec &Data, const ThreeDModelBase::tMeasPosVec &PosX,
        const ThreeDModelBase::tMeasPosVec &PosY,
        const ThreeDModelBase::tMeasPosVec &PosZ)
      {
        //do some consistency checks
        const size_t ndata = Data.size();
        assert(ndata == PosX.size());
        assert(ndata == PosY.size());
        assert(ndata == PosZ.size());
        std::ofstream outfile(filename.c_str());
        //first we have to write some general information about the file format
        WriteDataHeader(outfile, PosX, PosY, PosZ);

        outfile << "SCALARS " << DataName << " double" << std::endl;
        outfile << "LOOKUP_TABLE default" << std::endl;
        //and then just the data values
        std::copy(Data.begin(), Data.end(), std::ostream_iterator<double>(outfile, " "));
        outfile << std::endl;
        if (outfile.fail())
          throw FatalException("Problem writing vtk  file: " + filename);
      }

    /*! Write a collection of tensor measurements to a .vtk file for plotting
     * @param filename The name of the file, should contain the ending .vtk
     * @param DataName The name of the data for information in the plotting program
     * @param Data The vector of data tensors, 1 tensor per position
     * @param PosX The position of the measurement points in x-direction in m
     * @param PosY The position of the measurement points in y-direction in m
     * @param PosZ The position of the measurement points in z-direction in m
     */
    void Write3DTensorDataToVTK(const std::string &filename, const std::string &DataName,
        const jif3D::rvec &Data, const ThreeDModelBase::tMeasPosVec &PosX,
        const ThreeDModelBase::tMeasPosVec &PosY,
        const ThreeDModelBase::tMeasPosVec &PosZ)
      {
        //do some consistency checks
        const size_t ndata = Data.size();
        const size_t nmeas = PosX.size();
        assert(ndata == nmeas * 9);
        assert(ndata == PosY.size() * 9);
        assert(ndata == PosZ.size() * 9);

        std::ofstream outfile(filename.c_str());
        //first we have to write some general information about the file format
        WriteDataHeader(outfile, PosX, PosY, PosZ);

        outfile << "TENSORS " << DataName << " double" << std::endl;
        //and then just the data values
        //each line should contain three tensor components
        //and each tensor separated by two line breaks
        for (size_t i = 0; i < nmeas; ++i)
          {
            outfile << Data(i * 9) << " " << Data(i * 9 + 1) << " " << Data(i * 9 + 2)
                << "\n";
            outfile << Data(i * 9 + 3) << " " << Data(i * 9 + 4) << " " << Data(i * 9 + 5)
                << "\n";
            outfile << Data(i * 9 + 6) << " " << Data(i * 9 + 7) << " " << Data(i * 9 + 8)
                << "\n\n\n";
          }
        outfile << std::endl;
        if (outfile.fail())
          throw FatalException("Problem writing vtk  file: " + filename);
      }
  }
