/*
 * ThreeDMagnetizationModel.cpp
 *
 *  Created on: Apr 26, 2022
 *      Author: max
 */

#include "ThreeDMagnetizationModel.h"
#include "../ModelBase/VTKTools.h"
#include "../Global/NetCDFPortHelper.h"
#include "../ModelBase/NetCDFModelTools.h"
#include <netcdf>

using netCDF::NcFile;
using netCDF::NcVar;
using netCDF::NcDim;
using netCDF::NcVarAtt;

namespace jif3D
  {

    void ThreeDMagnetizationModel::WriteVTK(const std::string filename) const
      {
        Write3DVectorModelToVTK(filename, "Magnetization", GetXCoordinates(),
            GetYCoordinates(), GetZCoordinates(), GetData(), Magnetization_Y,
            Magnetization_Z);
      }
    void ThreeDMagnetizationModel::WriteNetCDF(const std::string &filename) const
      {
        netCDF::NcFile NetCDFFile(filename, netCDF::NcFile::replace);

        const ThreeDModelBase::t3DModelDim XCellCoords = GetXCoordinates();
        const ThreeDModelBase::t3DModelDim YCellCoords = GetYCoordinates();
        const ThreeDModelBase::t3DModelDim ZCellCoords = GetZCoordinates();


        //Make sure our data has the right dimensions and we have size information for each cell
        assert(GetData().num_dimensions() == 3);
        assert(GetData().shape()[0] == XCellCoords.size() - 1);
        assert(GetData().shape()[1] == YCellCoords.size() - 1);
        assert(GetData().shape()[2] == ZCellCoords.size() - 1);

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

        NcVar DataVar = NetCDFFile.addVar("Mx", netCDF::ncDouble, dimVec);
        DataVar.putAtt("units", "A/m");
        DataVar.putAtt("long_name", "X-component magnetization");
        NcVar VpDataVar = NetCDFFile.addVar("My", netCDF::ncDouble, dimVec);
        VpDataVar.putAtt("units", "A/m");
        VpDataVar.putAtt("long_name", "Y-component magnetization");
        NcVar DensDataVar = NetCDFFile.addVar("Mz", netCDF::ncDouble, dimVec);
        DensDataVar.putAtt("units", "A/m");
        DensDataVar.putAtt("long_name", "Z-component magnetization");
        //Write the model data itself
        //we store the data in the order Depth, Easting, Northing
        //but internally we work with x,y,z storage ordering
        //so we allocate a temporary buffer
        std::vector<double> databuffer(GetData().num_elements());
        std::vector<double> databuffer_y(Magnetization_Y.num_elements());
        std::vector<double> databuffer_z(Magnetization_Z.num_elements());
        //change the storage order
        const size_t xsize = XCoordDim.getSize();
        const size_t ysize = YCoordDim.getSize();
        const size_t zsize = ZCoordDim.getSize();

        for (size_t i = 0; i < zsize; ++i)
          for (size_t j = 0; j < ysize; ++j)
            for (size_t k = 0; k < xsize; ++k)
              {
                const size_t index = k + j * xsize + i * (xsize * ysize);
                databuffer[index] = GetData()[k][j][i];
                databuffer_y[index] = Magnetization_Y[k][j][i];
                databuffer_z[index] = Magnetization_Z[k][j][i];
              }
        //write the buffer to the file
        //        DataVar.put(databuffer, zsize, ysize, xsize);
        cxxport::put_legacy_ncvar(DataVar, databuffer.data(), zsize, ysize, xsize);
        cxxport::put_legacy_ncvar(VpDataVar, databuffer_y.data(), zsize, ysize, xsize);
        cxxport::put_legacy_ncvar(DensDataVar, databuffer_z.data(), zsize, ysize,
            xsize);

        NetCDFFile.putAtt("Conventions", "CF-1.3");
      }

    void ThreeDMagnetizationModel::ReadNetCDF(const std::string &model_file)
      {
        const std::string VsName = "Mx";
        const std::string VsUnit = "A/m";

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

        Magnetization_Y.resize(boost::extents[NX][NY][NZ]);
        Magnetization_Z.resize(boost::extents[NX][NY][NZ]);


        // Read model from nc file
        std::vector<double> tmp_data_y(NX * NY * NZ), tmp_data_z(NX * NY * NZ);
        NcVar densIn = ModelFile.getVar("My");
        densIn.getVar(tmp_data_y.data());
        NcVar vpIn = ModelFile.getVar("Mz");
        vpIn.getVar(tmp_data_z.data());

        for (int i = 0; i < NZ; ++i)
          {
            for (int j = 0; j < NY; ++j)
              {
                for (int k = 0; k < NX; ++k)
                  {
                    Magnetization_Y[k][j][i] = tmp_data_y[k + j * NX + i * (NX * NY)];
                    Magnetization_Z[k][j][i] = tmp_data_z[k + j * NX + i * (NX * NY)];
                  }
              }
          }
      }

    ThreeDMagnetizationModel::ThreeDMagnetizationModel()
      {
        // TODO Auto-generated constructor stub

      }

    ThreeDMagnetizationModel::~ThreeDMagnetizationModel()
      {
        // TODO Auto-generated destructor stub
      }

  } /* namespace jif3D */
