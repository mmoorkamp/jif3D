//============================================================================
// Name        : ThreeDModelBase.cpp
// Author      : Max Moorkamp
// Version     :
// Copyright   : 2008, MM
//============================================================================

#include "ThreeDModelBase.h"
#include "../Global/FatalException.h"
#include "NetCDFModelTools.h"
#include "VTKTools.h"
#include <cassert>
#include <fstream>
#include <boost/math/special_functions/relative_difference.hpp>

using netCDF::NcFile;
using netCDF::NcVar;
using netCDF::NcDim;

namespace jif3D
  {

    ThreeDModelBase::ThreeDModelBase() :
        Data(), XCellSizes(), YCellSizes(), ZCellSizes(), GridXCoordinates(), GridYCoordinates(), GridZCoordinates()
      {

      }

    ThreeDModelBase::ThreeDModelBase(const ThreeDModelBase& source) :
        Data(source.Data)
      {
        SetXCoordinates(source.GridXCoordinates);
        SetYCoordinates(source.GridYCoordinates);
        SetZCoordinates(source.GridZCoordinates);
      }

    ThreeDModelBase::~ThreeDModelBase()
      {

      }

    ThreeDModelBase& ThreeDModelBase::operator=(const ThreeDModelBase& source)
      {
        if (this == &source)
          return *this;
        //then we copy the data, i.e. the cells values of the 3D model
        Data.resize(
            boost::extents[source.Data.shape()[0]][source.Data.shape()[1]][source.Data.shape()[2]]);
        Data = source.Data;
        //now we copy the cell coordinates for all three directions
        //using the access function makes sure everything else is updated accordingly
        SetXCoordinates(source.GridXCoordinates);
        SetYCoordinates(source.GridYCoordinates);
        SetZCoordinates(source.GridZCoordinates);

        return *this;
      }
    bool ThreeDModelBase::operator ==(const ThreeDModelBase &b) const
      {
        double epsilon = 0.001;
        if (Data.num_elements() != b.Data.num_elements())
          return false;
        if (GridXCoordinates.size() != b.GridXCoordinates.size())
          return false;
        if (GridYCoordinates.size() != b.GridYCoordinates.size())
          return false;
        if (GridZCoordinates.size() != b.GridZCoordinates.size())
          return false;
        if (!std::equal(Data.origin(), Data.origin() + Data.num_elements(),
            b.Data.origin(), [epsilon](double a, double b)
              { return boost::math::relative_difference(a,b) < epsilon;}))
          return false;
        if (!std::equal(GridXCoordinates.begin(), GridXCoordinates.end(),
            b.GridXCoordinates.begin(), [epsilon](double a, double b)
              { return boost::math::relative_difference(a,b) < epsilon;}))
          return false;
        if (!std::equal(GridYCoordinates.begin(), GridYCoordinates.end(),
            b.GridYCoordinates.begin(), [epsilon](double a, double b)
              { return boost::math::relative_difference(a,b) < epsilon;}))
          return false;
        if (!std::equal(GridZCoordinates.begin(), GridZCoordinates.end(),
            b.GridZCoordinates.begin(), [epsilon](double a, double b)
              { return boost::math::relative_difference(a,b) < epsilon;}))
          return false;

        return true;
      }
    /*! Calculate the coordinates of the cells from the sizes and the Origin
     * @param Coordinates This vector will contain the coordinates of the left upper front corner of each cell
     * @param Sizes The size of each cell in m
     * @param ChangeFlag The flag that stores whether this coordinate has been changed
     * @param Origin The origin of the model axis
     */
    void ThreeDModelBase::CalcCoords(t3DModelDim &Coordinates, const t3DModelDim &Sizes)
      {
        if (!Sizes.empty())
          {
            //create a shorthand for the number of elements
            const size_t nelements = Sizes.size();
            //if the sizes have changed and there is something to calculate
            double Origin = Coordinates.size() >= 1 ? Coordinates.at(0) : 0;
            Coordinates.resize(nelements + 1);
            Coordinates[0] = Origin;
            double sum = Origin;
            for (size_t i = 0; i < nelements; ++i)
              {
                sum += Sizes[i];
                Coordinates[i + 1] = sum;
              }
          }
      }

    void ThreeDModelBase::CalcSizes(const t3DModelDim &Coordinates, t3DModelDim &Sizes)
      {
        if (!Coordinates.empty())
          {
            const size_t nelements = Coordinates.size();
            Sizes.resize(nelements - 1);
            for (size_t i = 1; i < nelements; ++i)
              {
                Sizes[i - 1] = Coordinates[i] - Coordinates[i - 1];
              }
          }
      }

    /*! For a given position within the model mesh find the three indices of the
     * corresponding cell in the mesh. This function is virtual as it provides
     * functionality for the general case where cells can have varying sizes.
     * For MT and seismic meshes we provide more efficient implementation.
     * @param xcoord The coordinate in x-direction in m
     * @param ycoord The coordinate in y-direction in m
     * @param zcoord The coordinate in z-direction in m
     * @return An array with the indices in x,y,z-direction, respectively
     */
    boost::array<ThreeDModelBase::t3DModelData::index, 3> ThreeDModelBase::FindAssociatedIndices(
        const double xcoord, const double ycoord, const double zcoord) const
      {
        //for all three directions we use the same approach
        //lower_bound gives an iterator to the field in Coordinates
        //that is just smaller then xcoord
        //the distance between this iterator and the beginning is the index
        const int xindex = std::distance(GetXCoordinates().begin(),
            std::lower_bound(GetXCoordinates().begin(), GetXCoordinates().end(), xcoord));
        const int yindex = std::distance(GetYCoordinates().begin(),
            std::lower_bound(GetYCoordinates().begin(), GetYCoordinates().end(), ycoord));
        const int zindex = std::distance(GetZCoordinates().begin(),
            std::lower_bound(GetZCoordinates().begin(), GetZCoordinates().end(), zcoord));
        //when we return the value we make sure that we cannot go out of bounds
        boost::array<t3DModelData::index, 3> idx =
              {
                { std::max(xindex - 1, 0), std::max(yindex - 1, 0), std::max(zindex - 1,
                    0) } };
        return idx;
      }

    void ThreeDModelBase::SetOrigin(const double x, const double y, const double z)
      {

        double OldX = GridXCoordinates.empty() ? 0.0 : GridXCoordinates[0];
        double OldY = GridYCoordinates.empty() ? 0.0 : GridYCoordinates[0];
        double OldZ = GridZCoordinates.empty() ? 0.0 : GridZCoordinates[0];
        //copy the information about the new origin
        std::transform(GridXCoordinates.begin(), GridXCoordinates.end(),
            GridXCoordinates.begin(), [x,OldX](double val)
              { return val +x - OldX;});
        std::transform(GridYCoordinates.begin(), GridYCoordinates.end(),
            GridYCoordinates.begin(), [y,OldY](double val)
              { return val +y - OldY;});
        std::transform(GridZCoordinates.begin(), GridZCoordinates.end(),
            GridZCoordinates.begin(), [z,OldZ](double val)
              { return val +z - OldZ;});
        //cell sizes do not change when shifting the origin, so nothing
        //needs to be done for those
      }

    void ThreeDModelBase::ReadDataFromNetCDF(const NcFile &NetCDFFile,
        const std::string &DataName, const std::string &UnitsName)
      {
        t3DModelDim XCD, YCD, ZCD;
        Read3DModelFromNetCDF(NetCDFFile, DataName, UnitsName, XCD, YCD, ZCD, Data);
        SetXCoordinates(XCD);
        SetYCoordinates(YCD);
        SetZCoordinates(ZCD);

        //we check that the sizes of the grid cell specifications and the data are matching
        if (GridXCoordinates.size() != Data.shape()[0] + 1
            || GridYCoordinates.size() != Data.shape()[1] + 1
            || GridZCoordinates.size() != Data.shape()[2] + 1)
          {
            throw jif3D::FatalException("Cell size specification does not match data !",
            __FILE__, __LINE__);
          }
      }

    void ThreeDModelBase::WriteDataToNetCDF(NcFile &NetCDFFile,
        const std::string &DataName, const std::string &UnitsName) const
      {
        Write3DModelToNetCDF(NetCDFFile, DataName, UnitsName, GetXCoordinates(),
            GetYCoordinates(), GetZCoordinates(), Data);
      }

    void ThreeDModelBase::WriteVTK(std::string filename,
        const std::string &DataName) const
      {
        Write3DModelToVTK(filename, DataName, GetXCoordinates(), GetYCoordinates(),
            GetZCoordinates(), GetData());
      }

    void ThreeDModelBase::WriteXYZ(const std::string &filename) const
      {
        std::ofstream outfile(filename.c_str());
        std::vector<double> XCenter(GridXCoordinates.size()), YCenter(
            GridYCoordinates.size()), ZCenter(GridZCoordinates.size());

        auto avgfunc = [](double a, double b)
          { return (a+b)/2.0;};
        std::adjacent_difference(GridXCoordinates.begin(), GridXCoordinates.end(),
            XCenter.begin(), avgfunc);
        std::adjacent_difference(GridYCoordinates.begin(), GridYCoordinates.end(),
            YCenter.begin(), avgfunc);
        std::adjacent_difference(GridZCoordinates.begin(), GridZCoordinates.end(),
            ZCenter.begin(), avgfunc);

        size_t ncells = GetNModelElements();
        for (size_t i = 0; i < ncells; ++i)
          {
            int xi, yi, zi;
            OffsetToIndex(i, xi, yi, zi);
            outfile << XCenter[xi + 1] << " " << YCenter[yi + 1] << " " << ZCenter[zi + 1]
                << " " << Data[xi][yi][zi] << "\n";
          }
      }

    void ThreeDModelBase::WriteUBC(const std::string &filename) const
      {
        std::ofstream meshfile((filename + ".mesh").c_str());
        meshfile << GetYCellSizes().size() << " ";
        meshfile << GetXCellSizes().size() << " ";
        meshfile << GetZCellSizes().size() << "\n";
        meshfile << GetYOrigin() << " " << GetXOrigin() << " " << GetZOrigin() << "\n";
        //so we have to swap x and y in the grid coordinate system
        std::copy(GetYCellSizes().begin(), GetYCellSizes().end(),
            std::ostream_iterator<double>(meshfile, " "));
        meshfile << "\n";
        std::copy(GetXCellSizes().begin(), GetXCellSizes().end(),
            std::ostream_iterator<double>(meshfile, " "));
        meshfile << "\n";
        std::copy(GetZCellSizes().begin(), GetZCellSizes().end(),
            std::ostream_iterator<double>(meshfile, " "));
        meshfile << "\n" << std::flush;

        std::ofstream valfile((filename + ".val").c_str());
        const size_t nx = GetXCellSizes().size();
        const size_t ny = GetYCellSizes().size();
        const size_t nz = GetZCellSizes().size();
        for (size_t i = 0; i < nx; ++i)
          {
            for (size_t j = 0; j < ny; ++j)
              {
                for (size_t k = 0; k < nz; ++k)
                  {
                    valfile << Data[i][j][k] << "\n";
                  }
              }
          }
        valfile << std::flush;

      }

  }
