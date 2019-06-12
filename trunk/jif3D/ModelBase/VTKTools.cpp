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
        const ThreeDModelBase::t3DModelDim &CellCoords)
      {
        //write the coordinate name, the number of cell boundary values and its type
        file << name << " " << CellCoords.size() << " double" << std::endl;
        std::copy(CellCoords.begin(),CellCoords.end(),std::ostream_iterator<double>(file," "));
        file << "\n";
      }

    //! A helper function to read the Coordinates for the 3 axes from a .vtk file
    void ReadCoordinatesFromVTK(std::ifstream &file,
        ThreeDModelBase::t3DModelDim &CellCoords)
      {
        //we need to read in a few key words from the vtk file
        //that are not of particular interest to us
        std::string dummy;
        int nvalues;
        //get the number of values and swallow up the vtk keywords
        file >> dummy >> nvalues >> dummy;
        CellCoords.resize(nvalues);
        //read in all the coordinate values
        for (double &val : CellCoords)
          {
            file >> val;
          }
      }

    /*! Write a .vtk file to plot a 3D model
     * @param filename The name of the output file, should contain the ending .vtk
     * @param DataName The name of the model data, for information for plotting programs
     * @param XCellCoords The coordinates of the cells in x-direction in m
     * @param YCellCoords The coordinates of the cells in y-direction in m
     * @param ZCellCoords The coordinates of the cells in z-direction in m
     * @param Data The model values within each cell, shape has to match the  cell sizes
     * @param XOrigin The x-origin of the coordinate system in m
     * @param YOrigin The y-origin of the coordinate system in m
     * @param ZOrigin The z-origin of the coordinate system in m
     */
    void Write3DModelToVTK(const std::string &filename, const std::string &DataName,
        const ThreeDModelBase::t3DModelDim &XCellCoords,
        const ThreeDModelBase::t3DModelDim &YCellCoords,
        const ThreeDModelBase::t3DModelDim &ZCellCoords,
        const ThreeDModelBase::t3DModelData &Data)
      {
        //do some consistency checks
        assert(Data.num_dimensions() == 3);
        assert(Data.shape()[0] == XCellCoords.size() - 1);
        assert(Data.shape()[1] == YCellCoords.size() - 1);
        assert(Data.shape()[2] == ZCellCoords.size() - 1);

        //get the size of each model direction
        const size_t nxvalues = XCellCoords.size() - 1;
        const size_t nyvalues = YCellCoords.size() - 1;
        const size_t nzvalues = ZCellCoords.size() - 1;

        std::ofstream outfile(filename.c_str());
        //first we have to write some general information about the file format
        outfile << "# vtk DataFile Version 2.0 \n";
        outfile << "3D Model data \n";
        outfile << "ASCII\n";
        outfile << "DATASET RECTILINEAR_GRID\n";
        //we write the left and right boundaries of each cell, but only store the left
        //so we have to write one extra value
        outfile << "DIMENSIONS " << nxvalues + 1 << " " << nyvalues + 1 << " "
            << nzvalues + 1 << "\n";
        //write information about the coordinate axes
        WriteCoordinatesToVTK(outfile, "X_COORDINATES", XCellCoords);
        WriteCoordinatesToVTK(outfile, "Y_COORDINATES", YCellCoords);
        WriteCoordinatesToVTK(outfile, "Z_COORDINATES", ZCellCoords);
        //write some information about the data itself
        outfile << "CELL_DATA " << nxvalues * nyvalues * nzvalues << "\n";
        outfile << "SCALARS " << DataName << " double" << "\n";
        outfile << "LOOKUP_TABLE default" << "\n";
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
          throw FatalException("Problem writing vtk  file: " + filename, __FILE__,
          __LINE__);
      }

    /*! Write a .vtk file to plot a 3D model with vector valued cells (.e.g the cross-gradient)
     * @param filename The name of the output file, should contain the ending .vtk
     * @param DataName The name of the model data, for information for plotting programs
     * @param XCellCoords The coordinates of the cells in x-direction in m
     * @param YCellCoords The coordinates of the cells in y-direction in m
     * @param ZCellCoords The coordinates of the cells in z-direction in m
     * @param CompX The x-component of the vector within each cell, shape has to match the  cell sizes
     * @param CompY The y-component of the vector within each cell, shape has to match the  cell sizes
     * @param CompZ The z-component of the vector within each cell, shape has to match the  cell sizes
     */
    void Write3DVectorModelToVTK(const std::string &filename, const std::string &DataName,
        const ThreeDModelBase::t3DModelDim &XCellCoords,
        const ThreeDModelBase::t3DModelDim &YCellCoords,
        const ThreeDModelBase::t3DModelDim &ZCellCoords,
        const ThreeDModelBase::t3DModelData &CompX,
        const ThreeDModelBase::t3DModelData &CompY,
        const ThreeDModelBase::t3DModelData &CompZ)
      {

        //get the size of each model direction
        const size_t nxvalues = XCellCoords.size() - 1;
        const size_t nyvalues = YCellCoords.size() - 1;
        const size_t nzvalues = ZCellCoords.size() - 1;

        std::ofstream outfile(filename.c_str());
        //first we have to write some general information about the file format
        outfile << "# vtk DataFile Version 2.0\n";
        outfile << "3D Model data\n";
        outfile << "ASCII\n";
        outfile << "DATASET RECTILINEAR_GRID\n";
        //we write the left and right boundaries of each cell, but only store the left
        //so we have to write one extra value
        outfile << "DIMENSIONS " << nxvalues + 1 << " " << nyvalues + 1 << " "
            << nzvalues + 1 << "\n";
        //write information about the coordinate axes
        WriteCoordinatesToVTK(outfile, "X_COORDINATES", XCellCoords);
        WriteCoordinatesToVTK(outfile, "Y_COORDINATES", YCellCoords);
        WriteCoordinatesToVTK(outfile, "Z_COORDINATES", ZCellCoords);
        //write some information about the data itself
        outfile << "CELL_DATA " << nxvalues * nyvalues * nzvalues << "\n";
        outfile << "VECTORS " << DataName << " double" << "\n";
        //and then just the data values
        for (size_t i = 0; i < nzvalues; ++i)
          {
            for (size_t j = 0; j < nyvalues; ++j)
              {
                for (size_t k = 0; k < nxvalues; ++k)
                  {
                    outfile << CompX[k][j][i] << " " << CompY[k][j][i] << " "
                        << CompZ[k][j][i] << "\n";
                  }
              }
          }
        if (outfile.fail())
          throw FatalException("Problem writing vtk  file: " + filename, __FILE__,
          __LINE__);
      }

    /*! Read a .vtk file containing a 3D model
     * @param filename The name of the input file, should contain the ending .vtk
     * @param XCellCoords The coordinates of the cells in x-direction in m
     * @param YCellCoords The coordinates of the cells in y-direction in m
     * @param ZCellCoords The coordinates of the cells in z-direction in m
     * @param Data The model values within each cell, shape has to match the  cell sizes
     */
    void Read3DModelFromVTK(const std::string &filename,
        ThreeDModelBase::t3DModelDim &XCellCoords,
        ThreeDModelBase::t3DModelDim &YCellCoords,
        ThreeDModelBase::t3DModelDim &ZCellCoords, ThreeDModelBase::t3DModelData &Data)
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
        ReadCoordinatesFromVTK(infile, XCellCoords);
        ReadCoordinatesFromVTK(infile, YCellCoords);
        ReadCoordinatesFromVTK(infile, ZCellCoords);
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
                      throw FatalException("Problem reading vtk  file: " + filename,
                      __FILE__, __LINE__);
                  }
              }
          }
      }


    /*! Read a .vtk file containing a 3D model
     * @param filename The name of the input file, should contain the ending .vtk
     * @param XCellCoords The coordinates of the cells in x-direction in m
     * @param YCellCoords The coordinates of the cells in y-direction in m
     * @param ZCellCoords The coordinates of the cells in z-direction in m
     * @param Data The model values within each cell, shape has to match the  cell sizes
     */
    void Read3DVectorModelFromVTK(const std::string &filename,
        ThreeDModelBase::t3DModelDim &XCellCoords,
        ThreeDModelBase::t3DModelDim &YCellCoords,
        ThreeDModelBase::t3DModelDim &ZCellCoords, ThreeDModelBase::t3DModelData &CompX,
         ThreeDModelBase::t3DModelData &CompY, ThreeDModelBase::t3DModelData &CompZ)
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
        ReadCoordinatesFromVTK(infile, XCellCoords);
        ReadCoordinatesFromVTK(infile, YCellCoords);
        ReadCoordinatesFromVTK(infile, ZCellCoords);
        //skip some more vtk specific information
        infile.getline(dummy, 255);
        infile.getline(dummy, 255);
        infile.getline(dummy, 255);
        CompX.resize(boost::extents[nxvalues][nyvalues][nzvalues]);
        CompY.resize(boost::extents[nxvalues][nyvalues][nzvalues]);
        CompZ.resize(boost::extents[nxvalues][nyvalues][nzvalues]);
        //read the values from the file, the storage ordering
        //is different from our internal storage, so we
        //have to use nested loops
        for (int i = 0; i < nzvalues; ++i)
          {
            for (int j = 0; j < nyvalues; ++j)
              {
                for (int k = 0; k < nxvalues; ++k)
                  {
                    infile >> CompX[k][j][i] >> CompY[k][j][i] >> CompZ[k][j][i];
                    if (infile.fail())
                      throw FatalException("Problem reading vtk  file: " + filename,
                      __FILE__, __LINE__);
                  }
              }
          }
      }

    //helper function that writes the common header information
    //for scalar and tensor data
    void WriteDataHeader(std::ofstream &outfile, const std::vector<double> &PosX,
        const std::vector<double> &PosY,
        const std::vector<double> &PosZ)
      {
        const size_t ndata = PosX.size();
        outfile << "# vtk DataFile Version 2.0\n";
        outfile << "3D Data\n";
        outfile << "ASCII\n";
        outfile << "DATASET UNSTRUCTURED_GRID\n";
        //we write the number of measurement points
        outfile << "POINTS " << ndata << " double \n";
        //write the coordinates of each point
        for (size_t i = 0; i < ndata; ++i)
          {
            outfile << PosX.at(i) << " " << PosY.at(i) << " " << PosZ.at(i) << "\n";
          }
        outfile << "POINT_DATA " << ndata << "\n";
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
        const std::vector<double> &Data, const std::vector<double> &PosX,
        const std::vector<double> &PosY,
        const std::vector<double> &PosZ)
      {
        //do some consistency checks
        const size_t ndata = Data.size();
        if (ndata != PosX.size())
          {
            throw FatalException(
                "Amount of x-coordinates does not match data size: " + filename, __FILE__,
                __LINE__);
          }
        if (ndata != PosY.size())
          {
            throw FatalException(
                "Amount of y-coordinates does not match data size: " + filename, __FILE__,
                __LINE__);
          }
        if (ndata != PosZ.size())
          {
            throw FatalException(
                "Amount of z-coordinates does not match data size: " + filename, __FILE__,
                __LINE__);
          }
        std::ofstream outfile(filename.c_str());
        //first we have to write some general information about the file format
        WriteDataHeader(outfile, PosX, PosY, PosZ);

        outfile << "SCALARS " << DataName << " double\n";
        outfile << "LOOKUP_TABLE default\n";
        //and then just the data values
        std::copy(Data.begin(), Data.end(), std::ostream_iterator<double>(outfile, " "));
        outfile << "\n";
        if (outfile.fail())
          throw FatalException("Problem writing vtk  file: " + filename, __FILE__,
          __LINE__);
      }

    /*! Write a collection of vector measurements to a .vtk file for plotting
     * @param filename The name of the file, should contain the ending .vtk
     * @param DataName The name of the data for information in the plotting program
     * @param Data The data, first the first component for all points, then the second then the third
     * @param PosX The position of the measurement points in x-direction in m
     * @param PosY The position of the measurement points in y-direction in m
     * @param PosZ The position of the measurement points in z-direction in m
     */
    void Write3DVectorDataToVTK(const std::string &filename, const std::string &DataName,
        const std::vector<double> &Data, const std::vector<double> &PosX,
        const std::vector<double> &PosY,
        const std::vector<double> &PosZ)
      {
        //do some consistency checks
        const size_t ndata = Data.size();
        const size_t nmeas = PosX.size();
        if (ndata != nmeas * 3)
          {
            throw FatalException(
                "Amount of x-coordinates does not match data size: " + filename, __FILE__,
                __LINE__);
          }
        if (nmeas != PosY.size())
          {
            throw FatalException(
                "Amount of y-coordinates does not match data size: " + filename, __FILE__,
                __LINE__);
          }
        if (nmeas != PosZ.size())
          {
            throw FatalException(
                "Amount of z-coordinates does not match data size: " + filename, __FILE__,
                __LINE__);
          }

        std::ofstream outfile(filename.c_str());
        //first we have to write some general information about the file format
        WriteDataHeader(outfile, PosX, PosY, PosZ);

        outfile << "VECTORS " << DataName << " double\n";
        //and then just the data values
        //each line should contain three tensor components
        //and each tensor separated by two line breaks
        for (size_t i = 0; i < nmeas; ++i)
          {
            outfile << Data.at(i * 9) << " " << Data.at(i * 9 + nmeas) << " "
                << Data.at(i * 9 + 2 * nmeas) << "\n";
          }
        outfile << "\n";
        if (outfile.fail())
          throw FatalException("Problem writing vtk  file: " + filename, __FILE__,
          __LINE__);
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
        const std::vector<double> &Data, const std::vector<double> &PosX,
        const std::vector<double> &PosY,
        const std::vector<double> &PosZ)
      {
        //do some consistency checks
        const size_t ndata = Data.size();
        const size_t nmeas = PosX.size();
        if (ndata != nmeas * 9)
          {
            throw FatalException(
                "Amount of x-coordinates does not match data size: " + filename, __FILE__,
                __LINE__);
          }
        if (nmeas != PosY.size())
          {
            throw FatalException(
                "Amount of y-coordinates does not match data size: " + filename, __FILE__,
                __LINE__);
          }
        if (nmeas != PosZ.size())
          {
            throw FatalException(
                "Amount of z-coordinates does not match data size: " + filename, __FILE__,
                __LINE__);
          }

        std::ofstream outfile(filename.c_str());
        //first we have to write some general information about the file format
        WriteDataHeader(outfile, PosX, PosY, PosZ);

        outfile << "TENSORS " << DataName << " double\n";
        //and then just the data values
        //each line should contain three tensor components
        //and each tensor separated by two line breaks
        for (size_t i = 0; i < nmeas; ++i)
          {
            outfile << Data.at(i * 9) << " " << Data.at(i * 9 + 1) << " " << Data.at(i * 9 + 2)
                << "\n";
            outfile << Data.at(i * 9 + 3) << " " << Data.at(i * 9 + 4) << " " << Data.at(i * 9 + 5)
                << "\n";
            outfile << Data.at(i * 9 + 6) << " " << Data.at(i * 9 + 7) << " " << Data.at(i * 9 + 8)
                << "\n\n\n";
          }
        outfile << "\n";
        if (outfile.fail())
          throw FatalException("Problem writing vtk  file: " + filename, __FILE__,
          __LINE__);
      }
  }
