/*
 * DCResistivityData.cpp
 *
 *  Created on: Nov 13, 2020
 *      Author: zhanjie
 */

#include "../Global/NetCDFTools.h"
#include "../Global/NetCDFPortHelper.h"
#include "../Global/VecMat.h"
#include "../ModelBase/NetCDFModelTools.h"
#include "../ModelBase/VTKTools.h"

#include <fstream>
#include <netcdf>
#include "DCResistivityData.h"

using netCDF::NcFile; // @suppress("Symbol is not resolved")
using netCDF::NcVar; // @suppress("Symbol is not resolved")
using netCDF::NcDim; // @suppress("Symbol is not resolved")
using netCDF::NcVarAtt; // @suppress("Symbol is not resolved")
using netCDF::ncInt; // @suppress("Symbol is not resolved")


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
        NcFile DataFile(filename, NcFile::read); // @suppress("Type cannot be resolved") // @suppress("Symbol is not resolved")
        //delete any old values in the model object
        ClearMeasurementPoints();
        ClearSecMeasurementPoint();
        ClearSourcePosPos();
        ClearSourceNegPos();
        ClearSourceIndices();

        //read the positions of the sources
        ReadVec(DataFile, DCSourcePosXName, SourcePosPosX); // @suppress("Invalid arguments")
        ReadVec(DataFile, DCSourcePosYName, SourcePosPosY); // @suppress("Invalid arguments")
        ReadVec(DataFile, DCSourcePosZName, SourcePosPosZ); // @suppress("Invalid arguments")
        ReadVec(DataFile, DCSourceNegXName, SourceNegPosX); // @suppress("Invalid arguments")
        ReadVec(DataFile, DCSourceNegYName, SourceNegPosY); // @suppress("Invalid arguments")
        ReadVec(DataFile, DCSourceNegZName, SourceNegPosZ); // @suppress("Invalid arguments")

        std::vector<double> MeasFirPosX, MeasFirPosY, MeasFirPosZ;
        //read the positions of the receivers
        ReadVec(DataFile, DCReceiverFirXName, MeasFirPosX); // @suppress("Invalid arguments")
        ReadVec(DataFile, DCReceiverFirYName, MeasFirPosY); // @suppress("Invalid arguments")
        ReadVec(DataFile, DCReceiverFirZName, MeasFirPosZ); // @suppress("Invalid arguments")
        ReadVec(DataFile, DCReceiverSecXName, MeasSecPosX); // @suppress("Invalid arguments")
        ReadVec(DataFile, DCReceiverSecYName, MeasSecPosY); // @suppress("Invalid arguments")
        ReadVec(DataFile, DCReceiverSecZName, MeasSecPosZ); // @suppress("Invalid arguments")


        //and configure the data object for these measurement point
        ReadVec(DataFile, DCSourceIndexName, SourceIndices); // @suppress("Invalid arguments")

        const size_t nmeas = MeasFirPosX.size();
//and add them as measurement positions to the data object
        for (size_t i = 0; i < nmeas; ++i)
          {
        	GeneralData::AddMeasurementPoint(MeasFirPosX[i], MeasFirPosY[i], MeasFirPosZ[i]);
          }

        //finally read in the apparent resistivities
        std::vector<double> Data, Error;
        ReadVec(DataFile, DCAppResistivityName, Data); // @suppress("Invalid arguments")

        // it is possible that there is no error information in the file
        // so we don't want to crash the program if not but set it to zero
        try
          {
            if (!DataFile.getVar(DCAppResistivityErrorName.c_str()).isNull()) // @suppress("Method cannot be resolved")
              {
                ReadVec(DataFile, DCAppResistivityErrorName, Error); // @suppress("Invalid arguments")
              }
            else
              {
                Error.resize(Data.size());
                std::fill(Error.begin(), Error.end(), 0.0);
              }
          } catch (netCDF::exceptions::NcException &ex) // @suppress("Type cannot be resolved")
          {
            // ignore
          }
        SetDataAndErrors(Data, Error);
      }


    void DCResistivityData::WriteNetCDF(const std::string &filename) const
      {
        //make sure all vectors have consistent sizes
        const size_t ndata = GetData().size();
        //assert(ndata == Model.GetMeasSecPosX().size());
        //ndata == Model.GetSourceIndices().size();

        const size_t nsourcepos = GetSourcePosPosX().size();
        //create a netcdf file
        NcFile DataFile(filename, NcFile::replace); // @suppress("Type cannot be resolved") // @suppress("Symbol is not resolved")

        //we use the source sequence number as a dimension
        NcDim SourceNumDim = DataFile.addDim(DCSourceNumberName, nsourcepos); // @suppress("Type cannot be resolved") // @suppress("Method cannot be resolved")

        //this is just an index over the source vector
        //and does not have any special meaning
        std::vector<int> SourcePosNumber(nsourcepos);
        std::iota(SourcePosNumber.begin(),SourcePosNumber.end(),0);
        NcVar SourceNumVar = DataFile.addVar(DCSourceNumberName, netCDF::ncInt, SourceNumDim); // @suppress("Type cannot be resolved") // @suppress("Symbol is not resolved") // @suppress("Method cannot be resolved")
        cxxport::put_legacy_ncvar(SourceNumVar, SourcePosNumber.data(), nsourcepos); // @suppress("Invalid arguments")

        //write out the source coordinates
        WriteVec(DataFile, DCSourcePosXName, GetSourcePosPosX(), SourceNumDim, "m"); // @suppress("Invalid arguments")
        WriteVec(DataFile, DCSourcePosYName, GetSourcePosPosY(), SourceNumDim, "m"); // @suppress("Invalid arguments")
        WriteVec(DataFile, DCSourcePosZName, GetSourcePosPosZ(), SourceNumDim, "m"); // @suppress("Invalid arguments")
        WriteVec(DataFile, DCSourceNegXName, GetSourceNegPosX(), SourceNumDim, "m"); // @suppress("Invalid arguments")
        WriteVec(DataFile, DCSourceNegYName, GetSourceNegPosY(), SourceNumDim, "m"); // @suppress("Invalid arguments")
        WriteVec(DataFile, DCSourceNegZName, GetSourceNegPosZ(), SourceNumDim, "m"); // @suppress("Invalid arguments")

        //write the index of the source for each measurement
        NcDim MeasIndexDim = DataFile.addDim(MeasIndexName, ndata); // @suppress("Type cannot be resolved") // @suppress("Method cannot be resolved")
        NcVar SourceIndexVar = DataFile.addVar(DCSourceIndexName, netCDF::ncInt, MeasIndexDim); // @suppress("Type cannot be resolved") // @suppress("Symbol is not resolved") // @suppress("Method cannot be resolved")
        cxxport::put_legacy_ncvar(SourceIndexVar, GetSourceIndices().data(), ndata); // @suppress("Invalid arguments")

        //write out the positions of the receivers, i.e. measurement positions
        WriteVec(DataFile, DCReceiverFirXName, GetMeasPosX(), MeasIndexDim, "m"); // @suppress("Invalid arguments")
        WriteVec(DataFile, DCReceiverFirYName, GetMeasPosY(), MeasIndexDim, "m"); // @suppress("Invalid arguments")
        WriteVec(DataFile, DCReceiverFirZName, GetMeasPosZ(), MeasIndexDim, "m"); // @suppress("Invalid arguments")
        WriteVec(DataFile, DCReceiverSecXName, GetMeasSecPosX(), MeasIndexDim, "m"); // @suppress("Invalid arguments")
        WriteVec(DataFile, DCReceiverSecYName, GetMeasSecPosY(), MeasIndexDim, "m"); // @suppress("Invalid arguments")
        WriteVec(DataFile, DCReceiverSecZName, GetMeasSecPosZ(), MeasIndexDim, "m"); // @suppress("Invalid arguments")

        //Write the apparent resistivity data and the error
        WriteVec(DataFile, DCAppResistivityName, GetData(), MeasIndexDim, "ohm.m"); // @suppress("Invalid arguments")
        WriteVec(DataFile, DCAppResistivityErrorName, GetErrors(), MeasIndexDim, "ohm.m"); // @suppress("Invalid arguments")

      }

    DCResistivityData::DCResistivityData()
      {

      }

    DCResistivityData::~DCResistivityData()
      {

      }


  } /* namespace jif3D */
