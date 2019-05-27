//============================================================================
// Name        : X3DModel.cpp
// Author      : Jul 2, 2009
// Version     :
// Copyright   : 2009, mmoorkamp
//============================================================================

#include "X3DModel.h"
#include "../Global/NumUtil.h"
#include <cassert>
#include <boost/numeric/conversion/cast.hpp>
#include <boost/math/special_functions/relative_difference.hpp>
#include <netcdf>
#include "../Global/NetCDFPortHelper.h"

using netCDF::NcFile;
using netCDF::NcVar;
using netCDF::NcDim;

namespace jif3D
  {

    X3DModel::X3DModel() :
        bg_thicknesses(), bg_conductivities()
      {

      }

    X3DModel::~X3DModel()
      {

      }



    X3DModel& X3DModel::operator=(const ThreeDModelBase& source)
      {
        if (&source != this)
          {
            ThreeDModelBase::operator=(source);
          }
        return *this;
      }

    bool X3DModel::operator ==(const X3DModel &b) const
      {
        double epsilon = 0.001;
        if (bg_thicknesses.size() != b.bg_thicknesses.size())
          return false;
        if (bg_conductivities.size() != b.bg_conductivities.size())
          return false;

        if (!std::equal(bg_thicknesses.begin(), bg_thicknesses.end(),
            b.bg_thicknesses.begin(), [epsilon](double a, double b)
              { return boost::math::relative_difference(a,b) < epsilon;}))
          return false;

        if (!std::equal(bg_conductivities.begin(), bg_conductivities.end(),
            b.bg_conductivities.begin(), [epsilon](double a, double b)
              { return boost::math::relative_difference(a,b) < epsilon;}))
          return false;

        //we only get here if nothing failed before
        return ThreeDModelBase::operator ==(b);
      }

    boost::array<ThreeDModelBase::t3DModelData::index, 3> X3DModel::FindAssociatedIndices(
        const double xcoord, const double ycoord, const double zcoord) const
      {
        const int xindex = boost::numeric_cast<int>(
            floor((xcoord - GetXCoordinates()[0]) / GetXCellSizes()[0]));
        const int yindex = boost::numeric_cast<int>(
            floor((ycoord - GetYCoordinates()[0]) / GetYCellSizes()[0]));
        const int zindex = std::distance(GetZCoordinates().begin(),
            std::lower_bound(GetZCoordinates().begin(), GetZCoordinates().end(), zcoord));
        //when we return the value we make sure that we cannot go out of bounds
        boost::array<t3DModelData::index, 3> idx =
          {
            { std::max(xindex, 0), std::max(yindex, 0), std::max(zindex - 1, 0) } };
        return idx;
      }

    void X3DModel::WriteNetCDF(const std::string &filename) const
      {
        NcFile DataFile(filename, NcFile::replace);
        //first write the 3D discretized part
        WriteDataToNetCDF(DataFile, ConductivityName, ConductivityUnit);
        //if we have some background layers, write them as well
        if (!bg_conductivities.empty())
          {
            assert(bg_conductivities.size() == bg_thicknesses.size());
            //we just number layers from 0 to n-1
            NcDim BackgroundDim = DataFile.addDim("bg_layers", bg_thicknesses.size());
            NcVar BackgroundVar = DataFile.addVar("bg_layers", netCDF::ncDouble,
                BackgroundDim);
            //here we generate the indices for the layers
            std::vector<double> layerindex(bg_thicknesses.size());
            std::iota(layerindex.begin(), layerindex.end(), 0);

//            BackgroundVar.put(&layerindex[0], layerindex.size());
            cxxport::put_legacy_ncvar(BackgroundVar, layerindex.data(),
                layerindex.size());
            BackgroundVar.putAtt("long_name", "Layer Index");
            //now we can write the actual parameters for the layers
            //first layer conductivities
            NcVar bgCondVar = DataFile.addVar("bg_conductivities", netCDF::ncDouble,
                BackgroundDim);
            bgCondVar.putAtt("long_name", "Background Conductivities");
            bgCondVar.putAtt("units", ConductivityUnit);
//            bgCondVar.put(&bg_conductivities[0], bg_conductivities.size());
            cxxport::put_legacy_ncvar(bgCondVar, bg_conductivities.data(),
                bg_conductivities.size());

            //and then layer thicknesses
            NcVar bgThickVar = DataFile.addVar("bg_thicknesses", netCDF::ncDouble,
                BackgroundDim);
            bgThickVar.putAtt("long_name", "Background Thicknesses");
            bgThickVar.putAtt("units", "m");

//            bgThickVar.put(&bg_thicknesses[0], bg_thicknesses.size());
            cxxport::put_legacy_ncvar(bgThickVar, bg_thicknesses.data(),
                bg_thicknesses.size());
          }
      }

    void X3DModel::ReadNetCDF(const std::string &filename)
      {
        //create the netcdf file object
        NcFile DataFile(filename, NcFile::read);
        //read in the 3D gridded data
        ReadDataFromNetCDF(DataFile, ConductivityName, ConductivityUnit);

        //change the error behaviour of the netcdf api
        try
          {
            //try to read the background data
            NcDim BackgroundDim = DataFile.getDim("bg_layers");
            //if we succeed
            if (!BackgroundDim.isNull())
              {
                //find out how many layers the background has
                const int nbglayers = BackgroundDim.getSize();
                //allocate memory
                bg_thicknesses.assign(nbglayers, 0.0);
                bg_conductivities.assign(nbglayers, 0.0);

                //and read in from the file
                NcVar bgDensVar = DataFile.getVar("bg_conductivities");
                NcVar bgThickVar = DataFile.getVar("bg_thicknesses");

//              bgDensVar.get(&bg_conductivities[0], nbglayers);
//              bgThickVar.get(&bg_thicknesses[0], nbglayers);

                cxxport::get_legacy_ncvar(bgDensVar, bg_conductivities.data(), nbglayers);
                cxxport::get_legacy_ncvar(bgThickVar, bg_thicknesses.data(), nbglayers);
              }
          } catch (const netCDF::exceptions::NcException &ex)
          {
            // ignore
          }
      }

  }
