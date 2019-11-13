/*
 * SurfaceWaveModels.cpp
 *
 *  Created on: 19 Sep 2019
 *      Author: bweise
 */

#include "../SurfaceWaves/SurfaceWaveModel.h"
#include "../Global/NetCDFPortHelper.h"
#include "../ModelBase/NetCDFModelTools.h"
#include <netcdf>
#include <vector>
using netCDF::NcFile;
using netCDF::NcVar;
using netCDF::NcDim;
using netCDF::NcVarAtt;

static const std::string OriginSuffix = "_Origin";

namespace jif3D
  {

    SurfaceWaveModel& SurfaceWaveModel::operator=(const ThreeDModelBase& source)
      {
        if (this == &source)
          return *this;
        ThreeDModelBase::operator =(source);

        return *this;
      }

    SurfaceWaveModel& SurfaceWaveModel::operator=(const SurfaceWaveModel& source)
      {
        if (this == &source)
          return *this;
        ThreeDModelBase::operator =(source);
        DataVp.resize(
            boost::extents[source.DataVp.shape()[0]][source.DataVp.shape()[1]][source.DataVp.shape()[2]]);
        DataVp = source.DataVp;
        DataDens.resize(
            boost::extents[source.DataDens.shape()[0]][source.DataDens.shape()[1]][source.DataDens.shape()[2]]);
        DataDens = source.DataDens;
        return *this;
      }

    void SurfaceWaveModel::ReadNetCDF(const std::string &model_file)
      {

        const std::string VsName = "Vs";
        const std::string VsUnit = "m/s";

        // Read density, vs, vp from nc file
        NcFile ModelFile(model_file, NcFile::read);

        int NX, NY, NZ;
        NcDim nxIn = ModelFile.getDim("Northing");
        NcDim nyIn = ModelFile.getDim("Easting");
        NcDim nzIn = ModelFile.getDim("Depth");
        NX = nxIn.getSize();
        NY = nyIn.getSize();
        NZ = nzIn.getSize();

        ReadDataFromNetCDF(ModelFile, VsName, VsUnit);

        DataVp.resize(boost::extents[NX][NY][NZ]);
        DataDens.resize(boost::extents[NX][NY][NZ]);

        // Read model from nc file
        std::vector<double> tmp_data_dns(NX * NY * NZ), tmp_data_vp(NX * NY * NZ);
        NcVar densIn = ModelFile.getVar("Density");
        densIn.getVar(tmp_data_dns.data());
        NcVar vpIn = ModelFile.getVar("Vp");
        vpIn.getVar(tmp_data_vp.data());

        for (int i = 0; i < NZ; ++i)
          {
            for (int j = 0; j < NY; ++j)
              {
                for (int k = 0; k < NX; ++k)
                  {
                    DataVp[k][j][i] = tmp_data_vp[k + j * NX + i * (NX * NY)];
                    DataDens[k][j][i] = tmp_data_dns[k + j * NX + i * (NX * NY)];
                  }
              }
          }
      }

    SurfaceWaveModel::SurfaceWaveModel()
      {
      }

    void SurfaceWaveModel::WriteNetCDF(const std::string &filename)
      {

        netCDF::NcFile NetCDFFile(filename, netCDF::NcFile::replace);

        const ThreeDModelBase::t3DModelDim XCellCoords = GetXCoordinates();
        const ThreeDModelBase::t3DModelDim YCellCoords = GetXCoordinates();
        const ThreeDModelBase::t3DModelDim ZCellCoords = GetXCoordinates();
        const ThreeDModelBase::t3DModelData Data = GetData();
        const ThreeDModelBase::t3DModelData VpData = GetVp();
        const ThreeDModelBase::t3DModelData DensData = GetDens();

        //Make sure our data has the right dimensions and we have size information for each cell
        assert(Data.num_dimensions() == 3);
        assert(Data.shape()[0] == XCellCoords.size() - 1);
        assert(Data.shape()[1] == YCellCoords.size() - 1);
        assert(Data.shape()[2] == ZCellCoords.size() - 1);

        //add information about boundary values
        NcDim BoundaryDim = NetCDFFile.addDim("nbound", 2);
        NcVar BoundaryVar = NetCDFFile.addVar("nbound", netCDF::ncInt, BoundaryDim);
        BoundaryVar.putAtt("long_name", "Boundary index variable");

        std::vector<int> BIndices;
        BIndices.push_back(1);
        BIndices.push_back(2);
        cxxport::put_legacy_ncvar(BoundaryVar, BIndices.data(), BIndices.size());
        //        BoundaryVar.put(&BIndices[0], 2); // old

        // Write the size information in x,y, and z-direction
        NcDim XCoordDim = WriteCoordinatesToNetCDF(NetCDFFile, "Northing", XCellCoords,
            BoundaryDim);
        NcDim YCoordDim = WriteCoordinatesToNetCDF(NetCDFFile, "Easting", YCellCoords,
            BoundaryDim);
        NcDim ZCoordDim = WriteCoordinatesToNetCDF(NetCDFFile, "Depth", ZCellCoords,
            BoundaryDim);

        std::vector<NcDim> dimVec;
        dimVec.push_back(ZCoordDim);
        dimVec.push_back(YCoordDim);
        dimVec.push_back(XCoordDim);

        NcVar DataVar = NetCDFFile.addVar("Vs", netCDF::ncDouble, dimVec);
        DataVar.putAtt("units", "m/s");
        DataVar.putAtt("long_name", "Shear wave velocity");
        NcVar VpDataVar = NetCDFFile.addVar("Vp", netCDF::ncDouble, dimVec);
        VpDataVar.putAtt("units", "m/s");
        VpDataVar.putAtt("long_name", "Pressure wave velocity");
        NcVar DensDataVar = NetCDFFile.addVar("Density", netCDF::ncDouble, dimVec);
        DensDataVar.putAtt("units", "kg/m3");
        DensDataVar.putAtt("long_name", "Density");
        //Write the model data itself
        //we store the data in the order Depth, Easting, Northing
        //but internally we work with x,y,z storage ordering
        //so we allocate a temporary buffer
        std::vector<double> databuffer(Data.num_elements());
        std::vector<double> databuffer_vp(VpData.num_elements());
        std::vector<double> databuffer_dens(DensData.num_elements());
        //change the storage order
        const size_t xsize = XCoordDim.getSize();
        const size_t ysize = YCoordDim.getSize();
        const size_t zsize = ZCoordDim.getSize();

        for (size_t i = 0; i < zsize; ++i)
          for (size_t j = 0; j < ysize; ++j)
            for (size_t k = 0; k < xsize; ++k)
              {
                databuffer[k + j * xsize + i * (xsize * ysize)] = Data[k][j][i];
                databuffer_vp[k + j * xsize + i * (xsize * ysize)] = VpData[k][j][i];
                databuffer_dens[k + j * xsize + i * (xsize * ysize)] = DensData[k][j][i];
              }
        //write the buffer to the file
        //        DataVar.put(databuffer, zsize, ysize, xsize);
        cxxport::put_legacy_ncvar(DataVar, databuffer.data(), zsize, ysize, xsize);
        cxxport::put_legacy_ncvar(VpDataVar, databuffer_vp.data(), zsize, ysize, xsize);
        cxxport::put_legacy_ncvar(DensDataVar, databuffer_dens.data(), zsize, ysize,
            xsize);

        NetCDFFile.putAtt("Conventions", "CF-1.3");
      }
  }

