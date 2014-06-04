//============================================================================
// Name        : ReadWriteDCResistivityData.h
// Author      : May 8, 2014
// Version     :
// Copyright   : 2014, zhanjie and mmoorkamp
//============================================================================

#include<fstream>
#include "ReadWriteDCResistivityData.h"
#include "../Global/NumUtil.h"
#include "../Global/NetCDFTools.h"
#include "../ModelBase/NetCDFModelTools.h"
#include "../ModelBase/ThreeDModelBase.h"

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

    void SaveApparentResistivity(const std::string &filename, const jif3D::rvec &Data,
        const jif3D::rvec &Error, const jif3D::ThreeDDCResistivityModel &Model)
      {
        //make sure all vectors have consistent sizes
        const size_t ndata = Data.size();
        //assert(ndata == Model.GetMeasSecPosX().size());
        //ndata == Model.GetSourceIndices().size();


        //create a netcdf file
        //create a netcdf file
        const size_t nsourcepos = Model.GetSourcePosPosX().size();
        const size_t nmeaspos = Model.GetMeasPosX().size();
        NcFile DataFile(filename.c_str(), NcFile::Replace);

        //we use the source sequence number as a dimension
        NcDim *SourceNumDim = DataFile.add_dim(DCSourceNumberName.c_str(), nsourcepos);

        //this is just an index over the source vector
        //and does not have any special meaning
        std::vector<int> SourcePosNumber;
        std::generate_n(std::back_inserter(SourcePosNumber), nsourcepos, IntSequence(0));
        NcVar *SourceNumVar = DataFile.add_var(DCSourceNumberName.c_str(), ncInt,
            SourceNumDim);
        SourceNumVar->put(&SourcePosNumber[0], nsourcepos);

        //write out the source coordinates
        WriteVec(DataFile, DCSourcePosXName, Model.GetSourcePosPosX(), SourceNumDim, "m");
        WriteVec(DataFile, DCSourcePosYName, Model.GetSourcePosPosY(), SourceNumDim, "m");
        WriteVec(DataFile, DCSourcePosZName, Model.GetSourcePosPosZ(), SourceNumDim, "m");
        WriteVec(DataFile, DCSourceNegXName, Model.GetSourceNegPosX(), SourceNumDim, "m");
        WriteVec(DataFile, DCSourceNegYName, Model.GetSourceNegPosY(), SourceNumDim, "m");
        WriteVec(DataFile, DCSourceNegZName, Model.GetSourceNegPosZ(), SourceNumDim, "m");

        //write the index of the source for each measurement
        NcDim *MeasIndexDim = DataFile.add_dim(MeasIndexName.c_str(), ndata);
        NcVar *SourceIndexVar = DataFile.add_var(DCSourceIndexName.c_str(), ncInt,
            MeasIndexDim);
        SourceIndexVar->put(&Model.GetSourceIndices()[0], ndata);

        //write out the positions of the receivers, i.e. measurement positions
        WriteVec(DataFile, DCReceiverFirXName, Model.GetMeasPosX(), MeasIndexDim, "m");
        WriteVec(DataFile, DCReceiverFirYName, Model.GetMeasPosY(), MeasIndexDim, "m");
        WriteVec(DataFile, DCReceiverFirZName, Model.GetMeasPosZ(), MeasIndexDim, "m");
        WriteVec(DataFile, DCReceiverSecXName, Model.GetMeasSecPosX(), MeasIndexDim, "m");
        WriteVec(DataFile, DCReceiverSecYName, Model.GetMeasSecPosY(), MeasIndexDim, "m");
        WriteVec(DataFile, DCReceiverSecZName, Model.GetMeasSecPosZ(), MeasIndexDim, "m");

        //Write the apparent resistivity data and the error
        WriteVec(DataFile, DCAppResistivityName, Data, MeasIndexDim, "ohm.m");
        WriteVec(DataFile, DCAppResistivityErrorName, Error, MeasIndexDim, "ohm.m");

      }

    void ReadApparentResistivity(const std::string &filename, jif3D::rvec &Data,
        jif3D::rvec &Error, jif3D::ThreeDDCResistivityModel &Model)
      {
        NcFile DataFile(filename.c_str(), NcFile::ReadOnly);
        jif3D::ThreeDModelBase::tMeasPosVec PosPosX, PosPosY, PosPosZ, NegPosX, NegPosY,
            NegPosZ, R1X, R1Y, R1Z, R2X, R2Y, R2Z;
//delete any old values in the model object
        Model.ClearMeasurementPoints();
        Model.ClearSecMeasurementPoint();
        Model.ClearSourcePosPos();
        Model.ClearSourceNegPos();
        Model.ClearSourceIndices();
//read the positions of the sources
        ReadVec(DataFile, DCSourcePosXName, PosPosX);
        ReadVec(DataFile, DCSourcePosYName, PosPosY);
        ReadVec(DataFile, DCSourcePosZName, PosPosZ);
        ReadVec(DataFile, DCSourceNegXName, NegPosX);
        ReadVec(DataFile, DCSourceNegYName, NegPosY);
        ReadVec(DataFile, DCSourceNegZName, NegPosZ);
        const size_t nsource = PosPosX.size();
//and add them to the model object
        for (size_t i = 0; i < nsource; ++i)
          {
            Model.AddSource(PosPosX[i], PosPosY[i], PosPosZ[i], NegPosX[i], NegPosY[i],
                NegPosZ[i]);
          }
//read the positions of the receivers
        ReadVec(DataFile, DCReceiverFirXName, R1X);
        ReadVec(DataFile, DCReceiverFirYName, R1Y);
        ReadVec(DataFile, DCReceiverFirZName, R1Z);
        ReadVec(DataFile, DCReceiverSecXName, R2X);
        ReadVec(DataFile, DCReceiverSecYName, R2Y);
        ReadVec(DataFile, DCReceiverSecZName, R2Z);

        std::vector<int> SourceIndices;
//now read the indices for the measurement position
        const size_t nconf = SourceIndices.size();
//and configure the model object for these measurement point
        ReadVec(DataFile, DCSourceIndexName, SourceIndices);

        const size_t nmeas = R1X.size();
//and add them as measurement positions to the model object
        for (size_t i = 0; i < nmeas; ++i)
          {
            Model.AddMeasurementPoint(R1X[i], R1Y[i], R1Z[i], R2X[i], R2Y[i], R2Z[i],
                SourceIndices[i]);
          }

//finally read in the traveltimes
        ReadVec(DataFile, DCAppResistivityName, Data);
// it is possible that there is no error information in the file
// so we don't want to crash the program if not but set it to zero
        NcError NetCDFError(NcError::silent_nonfatal);
        if (DataFile.get_var(DCAppResistivityErrorName.c_str()) != nullptr)
          {
            ReadVec(DataFile, DCAppResistivityErrorName, Error);
          }
        else
          {
            Error.resize(Data.size());
            Error.clear();
          }
      }

  }
