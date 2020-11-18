/*
 * DCResistivityData.cpp
 *
 *  Created on: Nov 13, 2020
 *      Author: zhanjie
 */

#include "../DCResistivity/DCResistivityData.h"

#include "../Global/NetCDFTools.h"
#include "../Global/NetCDFPortHelper.h"
#include "../Global/VecMat.h"
#include "../ModelBase/NetCDFModelTools.h"
#include "../ModelBase/VTKTools.h"

#include <fstream>
#include <netcdf>

using netCDF::NcFile;
using netCDF::NcVar;
using netCDF::NcDim;
using netCDF::NcVarAtt;
using netCDF::ncInt;


namespace jif3D
  {
    //we define these string variables to guarantee that we
    //use the exact same name for reading and writing the nectdf file
    static const std::string DCSourceNumberName = "DCSourceNumber";
    static const std::string DCSourcePosXName = "DCSourcePosPosX";
    static const std::string DCSourcePosYName = "DCSourcePosPosY";
    static const std::string DCSourcePosZName = "DCSourcePosPosZ";
    static const std::string DCSourceNegXName = "DCSourceNegPosX";
    static const std::string DCSourceNegYName = "DCSourceNegPosY";
    static const std::string DCSourceNegZName = "DCSourceNegPosZ";
    static const std::string DCReceiverFirXName = "DCReceiverFirX";
    static const std::string DCReceiverFirYName = "DCReceiverFirY";
    static const std::string DCReceiverFirZName = "DCReceiverFirZ";
    static const std::string DCReceiverSecXName = "DCReceiverSecX";
    static const std::string DCReceiverSecYName = "DCReceiverSecY";
    static const std::string DCReceiverSecZName = "DCReceiverSecZ";
    static const std::string DCSourceIndexName = "DCSourceIndex";
    static const std::string DCAppResistivityName = "DCAppResistivity";
    static const std::string DCAppResistivityErrorName = "DCErr";
    static const std::string MeasIndexName = "MeasIndex";

    void DCResistivityData::ReadNetCDF(const std::string &filename)
      {
        //create the netcdf file object
        NcFile DataFile(filename, NcFile::read);
        //delete any old values in the model object
        ClearMeasurementPoints();
        ClearSecMeasurementPoint();
        ClearSourcePosPos();
        ClearSourceNegPos();
        ClearSourceIndices();

        //read the positions of the sources
        ReadVec(DataFile, DCSourcePosXName, SourcePosPosX);
        ReadVec(DataFile, DCSourcePosYName, SourcePosPosY);
        ReadVec(DataFile, DCSourcePosZName, SourcePosPosZ);
        ReadVec(DataFile, DCSourceNegXName, SourceNegPosX);
        ReadVec(DataFile, DCSourceNegYName, SourceNegPosY);
        ReadVec(DataFile, DCSourceNegZName, SourceNegPosZ);
        const size_t nsource = SourcePosPosX.size();
//and add them to the data object
        for (size_t i = 0; i < nsource; ++i)
          {
            AddSource(SourcePosPosX[i], SourcePosPosY[i], SourcePosPosZ[i], SourceNegPosX[i], SourceNegPosY[i],
            		SourceNegPosZ[i]);
          }

        std::vector<double> MeasFirPosX, MeasFirPosY, MeasFirPosZ;
        //read the positions of the receivers
        ReadVec(DataFile, DCReceiverFirXName, MeasFirPosX);
        ReadVec(DataFile, DCReceiverFirYName, MeasFirPosY);
        ReadVec(DataFile, DCReceiverFirZName, MeasFirPosZ);
        ReadVec(DataFile, DCReceiverSecXName, MeasSecPosX);
        ReadVec(DataFile, DCReceiverSecYName, MeasSecPosY);
        ReadVec(DataFile, DCReceiverSecZName, MeasSecPosZ);


        //and configure the data object for these measurement point
        ReadVec(DataFile, DCSourceIndexName, SourceIndices);

        const size_t nmeas = MeasFirPosX.size();
//and add them as measurement positions to the data object
        for (size_t i = 0; i < nmeas; ++i)
          {
            AddMeasurementPoint(MeasFirPosX[i], MeasFirPosY[i], MeasFirPosZ[i], MeasSecPosX[i], MeasSecPosY[i], MeasSecPosZ[i],
                SourceIndices[i]);
          }

        //finally read in the apparent resistivities
        std::vector<double> Data, Error;
        ReadVec(DataFile, DCAppResistivityName, Data);

        // it is possible that there is no error information in the file
        // so we don't want to crash the program if not but set it to zero
        try
          {
            if (!DataFile.getVar(DCAppResistivityErrorName.c_str()).isNull())
              {
                ReadVec(DataFile, DCAppResistivityErrorName, Error);
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


    void DCResistivityData::WriteNetCDF(const std::string &filename)
      {
        //make sure all vectors have consistent sizes
        const size_t ndata = GetData().size();
        //assert(ndata == Model.GetMeasSecPosX().size());
        //ndata == Model.GetSourceIndices().size();

        const size_t nsourcepos = GetSourcePosPosX().size();
        //create a netcdf file
        NcFile DataFile(filename, NcFile::replace);

        //we use the source sequence number as a dimension
        NcDim SourceNumDim = DataFile.addDim(DCSourceNumberName, nsourcepos);

        //this is just an index over the source vector
        //and does not have any special meaning
        std::vector<int> SourcePosNumber(nsourcepos,0);
        std::iota(SourcePosNumber.begin(),SourcePosNumber.end(),0);

        NcVar SourceNumVar = DataFile.addVar(DCSourceNumberName, ncInt, SourceNumDim);
        cxxport::put_legacy_ncvar(SourceNumVar, SourcePosNumber.data(), nsourcepos);

        //write out the source coordinates
        WriteVec(DataFile, DCSourcePosXName, GetSourcePosPosX(), SourceNumDim, "m");
        WriteVec(DataFile, DCSourcePosYName, GetSourcePosPosY(), SourceNumDim, "m");
        WriteVec(DataFile, DCSourcePosZName, GetSourcePosPosZ(), SourceNumDim, "m");
        WriteVec(DataFile, DCSourceNegXName, GetSourceNegPosX(), SourceNumDim, "m");
        WriteVec(DataFile, DCSourceNegYName, GetSourceNegPosY(), SourceNumDim, "m");
        WriteVec(DataFile, DCSourceNegZName, GetSourceNegPosZ(), SourceNumDim, "m");

        //write the index of the source for each measurement
        NcDim MeasIndexDim = DataFile.addDim(MeasIndexName, ndata);
        NcVar SourceIndexVar = DataFile.addVar(DCSourceIndexName, ncInt, MeasIndexDim);
        cxxport::put_legacy_ncvar(SourceIndexVar, GetSourceIndices().data(), ndata);

        //write out the positions of the receivers, i.e. measurement positions
        WriteVec(DataFile, DCReceiverFirXName, GetMeasPosX(), MeasIndexDim, "m");
        WriteVec(DataFile, DCReceiverFirYName, GetMeasPosY(), MeasIndexDim, "m");
        WriteVec(DataFile, DCReceiverFirZName, GetMeasPosZ(), MeasIndexDim, "m");
        WriteVec(DataFile, DCReceiverSecXName, GetMeasSecPosX(), MeasIndexDim, "m");
        WriteVec(DataFile, DCReceiverSecYName, GetMeasSecPosY(), MeasIndexDim, "m");
        WriteVec(DataFile, DCReceiverSecZName, GetMeasSecPosZ(), MeasIndexDim, "m");

        //Write the apparent resistivity data and the error
        WriteVec(DataFile, DCAppResistivityName, GetData(), MeasIndexDim, "ohm.m");
        WriteVec(DataFile, DCAppResistivityErrorName, GetErrors(), MeasIndexDim, "ohm.m");

      }

    DCResistivityData::DCResistivityData()
      {

      }

    DCResistivityData::~DCResistivityData()
      {

      }


  } /* namespace jif3D */
