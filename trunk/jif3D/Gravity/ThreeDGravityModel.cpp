//============================================================================
// Name        : ThreeDGravityModel.cpp
// Author      : Max Moorkamp
// Version     :
// Copyright   : 2008, MM
//============================================================================

#include <cassert>

#include <fstream>
#include <iomanip>
#include <iterator>

#include <netcdf>

#include "ThreeDGravityModel.h"
#include "ReadWriteGravityData.h"
#include "../Global/NumUtil.h"
#include "../Global/NetCDFPortHelper.h"

using netCDF::NcFile;
using netCDF::NcVar;
using netCDF::NcDim;

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
        ThreeDModelBase(source), bg_densities(source.bg_densities), bg_thicknesses(
            source.bg_thicknesses)
      {

      }

    ThreeDGravityModel& ThreeDGravityModel::operator=(const ThreeDGravityModel& source)
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
            std::copy(source.bg_thicknesses.begin(), source.bg_thicknesses.end(),
                bg_thicknesses.begin());
          }
        return *this;
      }

    ThreeDGravityModel& ThreeDGravityModel::operator=(const ThreeDModelBase& source)
      {
        if (&source != this)
          {
            ThreeDModelBase::operator=(source);
          }
        return *this;
      }

    void ThreeDGravityModel::WriteNetCDF(const std::string &filename) const
      {
        NcFile DataFile(filename, NcFile::replace);
        //first write the 3D discretized part
        WriteDataToNetCDF(DataFile, DensityName, DensityUnit);
        //if we have some background layers, write them as well
        if (!bg_densities.empty())
          {
            assert(bg_densities.size() == bg_thicknesses.size());
            //Create the matching dimension for the background layers
            NcDim BackgroundDim = DataFile.addDim("bg_layers", bg_thicknesses.size());
            //now we can write the actual parameters for the layers
            NcVar bgDensVar = DataFile.addVar("bg_densities", netCDF::ncDouble, BackgroundDim);
            bgDensVar.putAtt("long_name", "Background Densities");
            bgDensVar.putAtt("units", DensityUnit);
            jif3D::cxxport::put_legacy_ncvar(bgDensVar, bg_densities.data(), bg_densities.size());
//            bgDensVar.put(&bg_densities[0], bg_densities.size()); // old

            NcVar bgThickVar = DataFile.addVar("bg_thicknesses", netCDF::ncDouble,
                BackgroundDim);
            bgThickVar.putAtt("long_name", "Background Thicknesses");
            bgThickVar.putAtt("units", "m");
            jif3D::cxxport::put_legacy_ncvar(bgThickVar, bg_thicknesses.data(), bg_thicknesses.size());
//            bgThickVar.put(&bg_thicknesses[0], bg_thicknesses.size()); // old
          }
      }

    void ThreeDGravityModel::ReadNetCDF(const std::string &filename)
      {

        //create the netcdf file object
        NcFile DataFile(filename, NcFile::read);
        //read in the 3D gridded data
        ReadDataFromNetCDF(DataFile, DensityName, DensityUnit);

        try {
          //try to read the background data
          NcDim BackgroundDim = DataFile.getDim("bg_layers");
          //if we succeed
          if (!BackgroundDim.isNull())
            {
              //find out how many layers the background has
              const size_t nbglayers = BackgroundDim.getSize();
              //allocate memory
              bg_densities.assign(nbglayers, 0.0);
              bg_thicknesses.assign(nbglayers, 0.0);
              //and read in from the file
              NcVar bgDensVar = DataFile.getVar("bg_densities");
              NcVar bgThickVar = DataFile.getVar("bg_thicknesses");

              jif3D::cxxport::get_legacy_ncvar(bgDensVar, bg_densities.data(), nbglayers);
              jif3D::cxxport::get_legacy_ncvar(bgThickVar, bg_thicknesses.data(), nbglayers);

//              bgDensVar.getVar(&bg_densities[0], nbglayers);
//              bgThickVar.getVar(&bg_thicknesses[0], nbglayers);
            }
        } catch(const netCDF::exceptions::NcException &ex) {
        	/*
        	 * Ignore exception. Before porting, the legacy NcError class was used to temporarily
        	 * set error handling to "silent_nonfatal".
        	 */
        }
      }

    void ThreeDGravityModel::ReadIgmas(const std::string &filename)
      {
        //the file we use to read the ascii data
        std::ifstream infile(filename.c_str());
        //first we read all values without any formatting into the vector values
        std::vector<double> values;
        std::copy(std::istream_iterator<double>(infile), std::istream_iterator<double>(),
            back_inserter(values));
        //if we don't have a multiple of 4 values there is a problem with the file
        if (values.size() % 4 != 0)
          {
            throw std::runtime_error("Problem reading file: " + filename);
          }
        const size_t nvalues = values.size() / 4;
        //some temporary object to store coordinates and density
        std::vector<double> xcoord(nvalues), ycoord(nvalues), zcoord(nvalues), density(
            nvalues);
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
        double minx = *std::min_element(xcoord.begin(), xcoord.end());
        std::transform(xcoord.begin(), xcoord.end(), xcoord.begin(), [minx] (double val)
          { return val - minx;});
        double miny = *std::min_element(ycoord.begin(), ycoord.end());
        std::transform(ycoord.begin(), ycoord.end(), ycoord.begin(), [miny] (double val)
          { return val - miny;});
        const size_t nx = xcoord.size();
        const size_t ny = ycoord.size();
        const size_t nz = zcoord.size();
        //check that the coordinate system and the densities are consistent
        if (nx * ny * nz != density.size())
          throw std::runtime_error("Mesh and density values are not consistent !");
        //reserve space and copy
        SetXCellSizes().resize(nx);
        SetYCellSizes().resize(ny);
        SetZCellSizes().resize(nz);
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
            SetDensities()[i % nx][(i / nx) % ny][nz - 1 - ((i / (nx * ny)) % nz)] =
                density.at(i);
          }
      }

  }
