//============================================================================
// Name        : X3DModel.cpp
// Author      : Jul 2, 2009
// Version     : 
// Copyright   : 2009, mmoorkamp
//============================================================================


#include "X3DModel.h"
#include "../Global/NumUtil.h"
namespace jiba
  {
    //we use these names for storing the model in NetCDF and vtk files
    static const std::string ConductivityName = "Conductivity";
    static const std::string ConductivityUnit = "S/m";

    X3DModel::X3DModel() :
      bg_thicknesses(), bg_conductivities()
      {

      }

    X3DModel::~X3DModel()
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
            std::copy(source.bg_thicknesses.begin(),
                source.bg_thicknesses.end(), bg_thicknesses.begin());
            std::copy(source.bg_conductivities.begin(),
                source.bg_conductivities.end(), bg_conductivities.begin());
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

    void X3DModel::WriteNetCDF(const std::string filename) const
      {
        NcFile DataFile(filename.c_str(), NcFile::Replace);
        //first write the 3D discretized part
        WriteDataToNetCDF(DataFile, ConductivityName, ConductivityUnit);
        //if we have some background layers, write them as well
        if (!bg_conductivities.empty())
          {
            assert(bg_conductivities.size() == bg_thicknesses.size());
            //we just number layers from 0 to n-1
            NcDim *BackgroundDim = DataFile.add_dim("bg_layers",
                bg_thicknesses.size());
            NcVar *BackgroundVar = DataFile.add_var("bg_layers", ncDouble,
                BackgroundDim);
            //here we generate the indices for the layers
            std::vector<double> layerindex;
            std::generate_n(back_inserter(layerindex), bg_thicknesses.size(),
                IntSequence(0));
            BackgroundVar->put(&layerindex[0], layerindex.size());
            BackgroundVar->add_att("long_name", "Layer Index");
            //now we can write the actual parameters for the layers
            //first layer conductivities
            NcVar *bgCondVar = DataFile.add_var("bg_conductivities", ncDouble,
                BackgroundDim);
            bgCondVar->add_att("long_name", "Background Conductivities");
            bgCondVar->add_att("units", ConductivityUnit.c_str());
            bgCondVar->put(&bg_conductivities[0], bg_conductivities.size());
            //and then layer thicknesses
            NcVar *bgThickVar = DataFile.add_var("bg_thicknesses", ncDouble,
                BackgroundDim);
            bgThickVar->add_att("long_name", "Background Thicknesses");
            bgThickVar->add_att("units", "m");
            bgThickVar->put(&bg_thicknesses[0], bg_thicknesses.size());
          }
      }

    void X3DModel::ReadNetCDF(const std::string filename)
      {
        //create the netcdf file object
        NcFile DataFile(filename.c_str(), NcFile::ReadOnly);
        //read in the 3D gridded data
        ReadDataFromNetCDF(DataFile, ConductivityName, ConductivityUnit);

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
            bg_thicknesses.assign(nbglayers, 0.0);
            bg_conductivities.assign(nbglayers, 0.0);
            //and read in from the file
            NcVar *bgDensVar = DataFile.get_var("bg_conductivities");
            NcVar *bgThickVar = DataFile.get_var("bg_thicknesses");
            bgDensVar->get(&bg_conductivities[0], nbglayers);
            bgThickVar->get(&bg_thicknesses[0], nbglayers);
          }
      }

  }
