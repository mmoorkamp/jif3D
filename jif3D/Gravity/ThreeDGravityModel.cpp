//============================================================================
// Name        : ThreeDGravityModel.cpp
// Author      : Max Moorkamp
// Version     :
// Copyright   : 2008, MM
//============================================================================


#include <boost/bind.hpp>
#include <cassert>
#include <netcdfcpp.h>
#include <fstream>
#include <iomanip>
#include <iterator>
#include "../Global/NumUtil.h"
#include "ThreeDGravityModel.h"
#include "ReadWriteGravityData.h"

namespace jif3D
  {
    static const std::string DensityName = "Density";
    static const std::string DensityUnit = "g/cm3";

    /*! The constructor does not take any parameters.
     */
    ThreeDGravityModel::ThreeDGravityModel() :
      bg_densities(), bg_thicknesses()
      {
      }

    ThreeDGravityModel::~ThreeDGravityModel()
      {
      }

    ThreeDGravityModel::ThreeDGravityModel(const ThreeDGravityModel &source) :
      ThreeDModelBase(source), bg_densities(source.bg_densities),
          bg_thicknesses(source.bg_thicknesses)
      {

      }

    ThreeDGravityModel& ThreeDGravityModel::operator=(
        const ThreeDGravityModel& source)
      {
        if (&source != this)
          {
            //in addition to copying the contents of the base class
            ThreeDModelBase::operator=(source);
            //we have to make sure we copy the information about the background
            bg_densities.resize(source.bg_densities.size());
            std::copy(source.bg_densities.begin(), source.bg_densities.end(),
                bg_densities.begin());
            bg_thicknesses.resize(source.bg_thicknesses.size());
            std::copy(source.bg_thicknesses.begin(),
                source.bg_thicknesses.end(), bg_thicknesses.begin());
          }
        return *this;
      }

    ThreeDGravityModel& ThreeDGravityModel::operator=(
        const ThreeDModelBase& source)
      {
        if (&source != this)
          {
            ThreeDModelBase::operator=(source);
          }
        return *this;
      }

    void ThreeDGravityModel::WriteNetCDF(const std::string filename) const
      {
        NcFile DataFile(filename.c_str(), NcFile::Replace);
        //first write the 3D discretized part
        WriteDataToNetCDF(DataFile, DensityName, DensityUnit);
        //if we have some background layers, write them as well
        if (!bg_densities.empty())
          {
            assert(bg_densities.size() == bg_thicknesses.size());
            //we just number layers from 0 to n-1
            NcDim *BackgroundDim = DataFile.add_dim("bg_layers",
                bg_thicknesses.size());
            NcVar *BackgroundVar = DataFile.add_var("bg_layers", ncDouble,
                BackgroundDim);
            std::vector<double> layerindex;
            std::generate_n(back_inserter(layerindex), bg_thicknesses.size(),
                IntSequence(0));
            BackgroundVar->put(&layerindex[0], layerindex.size());
            BackgroundVar->add_att("long_name", "Layer Index");
            //now we can write the actual parameters for the layers
            NcVar *bgDensVar = DataFile.add_var("bg_densities", ncDouble,
                BackgroundDim);
            bgDensVar->add_att("long_name", "Background Densities");
            bgDensVar->add_att("units", DensityUnit.c_str());
            NcVar *bgThickVar = DataFile.add_var("bg_thicknesses", ncDouble,
                BackgroundDim);
            bgThickVar->add_att("long_name", "Background Thicknesses");
            bgThickVar->add_att("units", "m");
            bgDensVar->put(&bg_densities[0], bg_densities.size());
            bgThickVar->put(&bg_thicknesses[0], bg_thicknesses.size());
          }
      }

    void ThreeDGravityModel::ReadNetCDF(const std::string filename)
      {

        //create the netcdf file object
        NcFile DataFile(filename.c_str(), NcFile::ReadOnly);
        //read in the 3D gridded data
        ReadDataFromNetCDF(DataFile, DensityName, DensityUnit);

        //change the error behaviour of the netcdf api
        NcError err(NcError::silent_nonfatal);
        //try to read the background data
        NcDim *BackgroundDim = DataFile.get_dim("bg_layers");
        //if we succeed
        if (BackgroundDim)
          {
            //find out how many layers the background has
            const int nbglayers = BackgroundDim->size();
            //allocate memory
            bg_densities.assign(nbglayers, 0.0);
            bg_thicknesses.assign(nbglayers, 0.0);
            //and read in from the file
            NcVar *bgDensVar = DataFile.get_var("bg_densities");
            NcVar *bgThickVar = DataFile.get_var("bg_thicknesses");
            bgDensVar->get(&bg_densities[0], nbglayers);
            bgThickVar->get(&bg_thicknesses[0], nbglayers);
          }
      }

    void ThreeDGravityModel::ReadIgmas(const std::string filename)
      {
        //the file we use to read the ascii data
        std::ifstream infile(filename.c_str());
        //first we read all values without any formatting into the vector values
        std::vector<double> values;
        std::copy(std::istream_iterator<double>(infile), std::istream_iterator<
            double>(), back_inserter(values));
        //if we don't have a multiple of 4 values there is a problem with the file
        if (values.size() % 4 != 0)
          {
            throw std::runtime_error("Problem reading file: " + filename);
          }
        const size_t nvalues = values.size() / 4;
        //some temporary object to store coordinates and density
        std::vector<double> xcoord(nvalues), ycoord(nvalues), zcoord(nvalues),
            density(nvalues);
        //first we copy without consideration of storage order and duplicates
        for (size_t i = 0; i < values.size(); i += 4)
          {
            const size_t index = i / 4;
            xcoord.at(index) = values.at(i);
            ycoord.at(index) = values.at(i + 1);
            zcoord.at(index) = values.at(i + 2);
            density.at(index) = values.at(i + 3);
          }
        //now we sort out :-) the actual coordinates
        std::sort(xcoord.begin(), xcoord.end());
        std::sort(ycoord.begin(), ycoord.end());
        std::sort(zcoord.begin(), zcoord.end());
        xcoord.erase(std::unique(xcoord.begin(), xcoord.end()), xcoord.end());
        ycoord.erase(std::unique(ycoord.begin(), ycoord.end()), ycoord.end());
        zcoord.erase(std::unique(zcoord.begin(), zcoord.end()), zcoord.end());
        //our z-coordinates are positive down, here they are positive up
        std::transform(zcoord.begin(), zcoord.end(), zcoord.begin(),
            std::negate<double>());
        //we always start with 0, in igmas it can be any number
        std::transform(xcoord.begin(), xcoord.end(), xcoord.begin(),
            boost::bind(std::minus<double>(), _1, *std::min_element(
                xcoord.begin(), xcoord.end())));
        std::transform(ycoord.begin(), ycoord.end(), ycoord.begin(),
            boost::bind(std::minus<double>(), _1, *std::min_element(
                ycoord.begin(), ycoord.end())));
        const size_t nx = xcoord.size();
        const size_t ny = ycoord.size();
        const size_t nz = zcoord.size();
        //check that the coordinate system and the densities are consistent
        if (nx * ny * nz != density.size())
          throw std::runtime_error(
              "Mesh and density values are not consistent !");
        //reserve space and copy
        SetXCellSizes().resize(boost::extents[nx]);
        SetYCellSizes().resize(boost::extents[ny]);
        SetZCellSizes().resize(boost::extents[nz]);
        SetDensities().resize(boost::extents[nx][ny][nz]);
        //we skip the first one, because it is always zero
        std::adjacent_difference(xcoord.begin() + 1, xcoord.end(),
            SetXCellSizes().begin());
        std::adjacent_difference(ycoord.begin() + 1, ycoord.end(),
            SetYCellSizes().begin());
        std::adjacent_difference(zcoord.rbegin() + 1, zcoord.rend(),
            SetZCellSizes().begin());
        //we need z varying fastest, but it is x varying fastest
        for (size_t i = 0; i < nvalues; ++i)
          {
            SetDensities()[i % nx][(i / nx) % ny][nz - 1 - ((i / (nx * ny))
                % nz)] = density.at(i);
          }
      }

  }
