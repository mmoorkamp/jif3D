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

    X3DModel::X3DModel(const X3DModel &source) :
        ThreeDMTModel(source), bg_thicknesses(source.bg_thicknesses), bg_conductivities(
            source.bg_conductivities)
      {

      }

    X3DModel& X3DModel::operator=(const X3DModel& source)
      {
        if (&source != this)
          {
            //first we copy the base class
            ThreeDMTModel::operator=(source);
            //then we copy the additional information about the background layers
            //that is not contained in the base class
            bg_thicknesses.resize(source.bg_thicknesses.size());
            bg_conductivities.resize(source.bg_conductivities.size());
            std::copy(source.bg_thicknesses.begin(), source.bg_thicknesses.end(),
                bg_thicknesses.begin());
            std::copy(source.bg_conductivities.begin(), source.bg_conductivities.end(),
                bg_conductivities.begin());
          }
        return *this;
      }

    X3DModel& X3DModel::operator=(const ThreeDModelBase& source)
      {
        if (&source != this)
          {
            ThreeDModelBase::operator=(source);
          }
        return *this;
      }

    boost::array<ThreeDModelBase::t3DModelData::index, 3> X3DModel::FindAssociatedIndices(
        const double xcoord, const double ycoord, const double zcoord) const
      {
        const int xindex = boost::numeric_cast<int>(
            floor((xcoord - XOrigin) / GetXCellSizes()[0]));
        const int yindex = boost::numeric_cast<int>(
            floor((ycoord - YOrigin) / GetYCellSizes()[0]));
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
            NcVar BackgroundVar = DataFile.addVar("bg_layers", netCDF::ncDouble, BackgroundDim);
            //here we generate the indices for the layers
            std::vector<double> layerindex;
            std::generate_n(back_inserter(layerindex), bg_thicknesses.size(),
                IntSequence(0));

//            BackgroundVar.put(&layerindex[0], layerindex.size());
            cxxport::put_legacy_ncvar(BackgroundVar, layerindex.data(), layerindex.size());
            BackgroundVar.putAtt("long_name", "Layer Index");
            //now we can write the actual parameters for the layers
            //first layer conductivities
            NcVar bgCondVar = DataFile.addVar("bg_conductivities", netCDF::ncDouble,
                BackgroundDim);
            bgCondVar.putAtt("long_name", "Background Conductivities");
            bgCondVar.putAtt("units", ConductivityUnit);
//            bgCondVar.put(&bg_conductivities[0], bg_conductivities.size());
            cxxport::put_legacy_ncvar(bgCondVar, bg_conductivities.data(), bg_conductivities.size());

            //and then layer thicknesses
            NcVar bgThickVar = DataFile.addVar("bg_thicknesses", netCDF::ncDouble,
                BackgroundDim);
            bgThickVar.putAtt("long_name", "Background Thicknesses");
            bgThickVar.putAtt("units", "m");

//            bgThickVar.put(&bg_thicknesses[0], bg_thicknesses.size());
            cxxport::put_legacy_ncvar(bgThickVar, bg_thicknesses.data(), bg_thicknesses.size());
          }
      }

    void X3DModel::ReadNetCDF(const std::string &filename)
      {
        //create the netcdf file object
        NcFile DataFile(filename, NcFile::read);
        //read in the 3D gridded data
        ReadDataFromNetCDF(DataFile, ConductivityName, ConductivityUnit);

        //change the error behaviour of the netcdf api
        try {
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
        } catch(const netCDF::exceptions::NcException &ex) {
          // ignore
        }
      }

  }
