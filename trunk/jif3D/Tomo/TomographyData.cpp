/*
 * TomographyData.cpp
 *
 *  Created on: May 21, 2019
 *      Author: max
 */

#include "../Tomo/TomographyData.h"

#include "../Global/NetCDFTools.h"
#include "../Global/NetCDFPortHelper.h"
#include "../Global/VecMat.h"
#include "../ModelBase/NetCDFModelTools.h"
#include "../ModelBase/VTKTools.h"

#include <fstream>
#include <netcdf>

using netCDF::NcFile;
using netCDF::NcDim;
using netCDF::NcVar;
namespace jif3D
  {

    //we define these string variables to guarantee that we
    //use the exact same name for reading and writing the nectdf file
    static const std::string SourceNumberName = "SourceNumber";
    static const std::string ReceiverNumberName = "ReceiverNumber";
    static const std::string SourcePosXName = "SourcePosX";
    static const std::string SourcePosYName = "SourcePosY";
    static const std::string SourcePosZName = "SourcePosZ";
    static const std::string MeasIndexName = "MeasIndex";
    static const std::string TravelTimeName = "TravelTime";
    static const std::string TravelTimeErrorName = "dT";
    static const std::string SourceIndexName = "SourceIndex";
    static const std::string ReceiverIndexName = "ReceiverIndex";

    void TomographyData::ReadNetCDF(const std::string &filename)
      {
        NcFile DataFile(filename, NcFile::read);
        //delete any old values in the model object
        ClearMeasurementPoints();
        ClearSourcePos();
        ClearMeasurementConfigurations();
        //read the positions of the sources
        ReadVec(DataFile, SourcePosXName, SourcePosX);
        ReadVec(DataFile, SourcePosYName, SourcePosY);
        ReadVec(DataFile, SourcePosZName, SourcePosZ);

        std::vector<double> PosX, PosY, PosZ;
        //read the positions of the receivers
        ReadVec(DataFile, MeasPosXName, PosX);
        ReadVec(DataFile, MeasPosYName, PosY);
        ReadVec(DataFile, MeasPosZName, PosZ);
        SetMeasurementPoints(PosX, PosY, PosZ);

        //now read the indices for the source receiver combinations
        //for each measurement
        ReadVec(DataFile, SourceIndexName, SourceIndices);
        ReadVec(DataFile, ReceiverIndexName, ReceiverIndices);
        //finally read in the traveltimes
        std::vector<double> Data, Error;
        ReadVec(DataFile, TravelTimeName, Data);
        // it is possible that there is no error information in the file
        // so we don't want to crash the program if not but set it to zero
        try
          {
            if (!DataFile.getVar(TravelTimeErrorName.c_str()).isNull())
              {
                ReadVec(DataFile, TravelTimeErrorName, Error);
              }
            else
              {
                Error.resize(Data.size());
                std::fill(Error.begin(), Error.end(), 0.0);
              }
          } catch (netCDF::exceptions::NcException &ex)
          {
            // ignore
          }
        SetDataAndErrors(Data, Error);
      }
    void TomographyData::WriteNetCDF(const std::string &filename)
      {
        //make sure all vectors have consistent sizes
        const size_t ndata = GetData().size();
        assert(ndata == GetReceiverIndices().size());
        assert(ndata == GetSourceIndices().size());

        //create a netcdf file
        const size_t nsourcepos = GetSourcePosX().size();
        const size_t nrecpos = GetMeasPosX().size();
        NcFile DataFile(filename, NcFile::replace);
        //we use the station number as a dimension
        NcDim SourceNumDim = DataFile.addDim(SourceNumberName, nsourcepos);
        NcDim RecNumDim = DataFile.addDim(ReceiverNumberName, nrecpos);

        //write out the measurement coordinates
        WriteVec(DataFile, SourcePosXName, GetSourcePosX(), SourceNumDim, "m");
        WriteVec(DataFile, SourcePosYName, GetSourcePosY(), SourceNumDim, "m");
        WriteVec(DataFile, SourcePosZName, GetSourcePosZ(), SourceNumDim, "m");
        //write out the positions of the receivers, i.e. measurement positions
        WriteVec(DataFile, MeasPosXName, GetMeasPosX(), RecNumDim, "m");
        WriteVec(DataFile, MeasPosYName, GetMeasPosY(), RecNumDim, "m");
        WriteVec(DataFile, MeasPosZName, GetMeasPosZ(), RecNumDim, "m");
        //generate an index for the receivers
        NcDim MeasIndexDim = DataFile.addDim(MeasIndexName, ndata);
        std::vector<int> MeasIndex(ndata);
        std::iota(MeasIndex.begin(), MeasIndex.end(), 0);
        NcVar MeasIndexVar = DataFile.addVar(MeasIndexName, netCDF::ncInt, MeasIndexDim);

        //        MeasIndexVar.put(&MeasIndex[0], ndata);
        cxxport::put_legacy_ncvar(MeasIndexVar, MeasIndex.data(), ndata);
        //Write the travel times and the error
        WriteVec(DataFile, TravelTimeName, GetData(), MeasIndexDim, "s");

        //write the index of the source for each measurement
        NcVar SourceIndexVar = DataFile.addVar(SourceIndexName, netCDF::ncInt,
            MeasIndexDim);
        //        SourceIndexVar.put(&Model.GetSourceIndices()[0], MeasIndexDim.getSize());
        cxxport::put_legacy_ncvar(SourceIndexVar, GetSourceIndices().data(),
            MeasIndexDim.getSize());
        //write the index of the receiver for each measurement
        NcVar RecIndexVar = DataFile.addVar(ReceiverIndexName, netCDF::ncInt,
            MeasIndexDim);

        //        RecIndexVar.put(&Model.GetReceiverIndices()[0], MeasIndexDim.getSize());
        cxxport::put_legacy_ncvar(RecIndexVar, GetReceiverIndices().data(),
            MeasIndexDim.getSize());

        WriteVec(DataFile, TravelTimeErrorName, GetErrors(), MeasIndexDim, "s");

      }

    void TomographyData::WriteSourcePoints(const std::string &filename)
      {
        std::vector<double> SourceNum(GetSourcePosX().size());
        std::iota(SourceNum.begin(), SourceNum.end(), 1);
        jif3D::Write3DDataToVTK(filename, "Sources", SourceNum, GetSourcePosX(),
            GetSourcePosY(), GetSourcePosZ());
      }

    TomographyData::TomographyData()
      {
        // TODO Auto-generated constructor stub

      }

    TomographyData::~TomographyData()
      {
        // TODO Auto-generated destructor stub
      }

  } /* namespace jif3D */
