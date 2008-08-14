//============================================================================
// Name        : ReadWriteGravityData.cpp
// Author      : Max Moorkamp
// Version     :
// Copyright   : 2008, MM
//============================================================================

#include "ReadWriteGravityData.h"
#include <netcdfcpp.h>
#include "../Global/NumUtil.h"
#include <boost/bind.hpp>

namespace jiba
  {
    static const std::string MeasPosXName = "MeasPosX";
    static const std::string MeasPosYName = "MeasPosY";
    static const std::string MeasPosZName = "MeasPosZ";
    static const std::string ScalarGravityName = "Scalar_gravity";
    static const std::string StationNumberName = "StationNumber";
    static const std::string UxxName = "Uxx";
    static const std::string UxyName = "Uxy";
    static const std::string UxzName = "Uxz";
    static const std::string UyxName = "Uyx";
    static const std::string UyyName = "Uyy";
    static const std::string UyzName = "Uyz";
    static const std::string UzxName = "Uzx";
    static const std::string UzyName = "Uzy";
    static const std::string UzzName = "Uzz";

    //! Read one measurement position coordinate from a netcdf file
    void ReadVec(NcFile &NetCDFFile, const std::string &MeasPosName,
        ThreeDGravityModel::tMeasPosVec &Position)
      {
        //create a netcdf dimension for the Station number
        NcDim *Dim = NetCDFFile.get_dim(StationNumberName.c_str());
        //determine the size of that dimension
        const size_t nvalues = Dim->size();

        //allocate memory in the class variable
        Position.assign(nvalues, 0.0);
        // create netcdf variable with the same name as the dimension
        NcVar *SizeVar = NetCDFFile.get_var(MeasPosName.c_str());
        //read coordinate values from netcdf file
        SizeVar->get(&Position[0], nvalues);
      }

    void ReadMatComp(NcFile &NetCDFFile, const std::string &CompName,
        ThreeDGravityModel::tTensorMeasVec &MatVec, const size_t n,
        const size_t m)
      {
        ThreeDGravityModel::tScalarMeasVec tempdata(MatVec.size());
        ReadVec(NetCDFFile, CompName, tempdata);
        assert(tempdata.size() == MatVec.size());
        for (size_t i = 0; i < tempdata.size(); ++i)
          MatVec.at(i)(n, m) = tempdata.at(i);
      }

    void WriteVec(NcFile &NetCDFFile, const std::string &MeasPosName,
        const ThreeDGravityModel::tMeasPosVec &Position, NcDim *Dimension)
      {
        const size_t nmeas = Position.size();
        NcVar *PosVar = NetCDFFile.add_var(MeasPosName.c_str(), ncDouble,
            Dimension);
        PosVar->add_att("units", "m");
        PosVar->put(&Position[0], nmeas);
      }

    void WriteMatComp(NcFile &NetCDFFile, const std::string &CompName,
        ThreeDGravityModel::tTensorMeasVec &MatVec, const size_t n,
        const size_t m, NcDim *Dimension)
      {
        ThreeDGravityModel::tScalarMeasVec tempdata(MatVec.size());
        for (size_t i = 0; i < tempdata.size(); ++i)
          tempdata.at(i) = MatVec.at(i)(n, m);
        WriteVec(NetCDFFile, CompName, tempdata, Dimension);
      }

    void MakeStatNumDim(const size_t nmeas, NcDim *Dimension,NcFile &DataFile )
      {
        std::vector<int> StationNumber;
        std::generate_n(back_inserter(StationNumber), nmeas, IntSequence(0));
        Dimension = DataFile.add_dim(StationNumberName.c_str(), nmeas);
        NcVar *StatNumVar = DataFile.add_var(StationNumberName.c_str(), ncInt,
            Dimension);
        StatNumVar->put(&StationNumber[0], nmeas);
      }

    void ReadMeasPosNetCDF(const std::string filename,
        ThreeDGravityModel::tMeasPosVec &PosX,
        ThreeDGravityModel::tMeasPosVec &PosY,
        ThreeDGravityModel::tMeasPosVec &PosZ)
      {
        NcFile DataFile(filename.c_str(), NcFile::ReadOnly);
        ReadVec(DataFile, MeasPosXName, PosX);
        ReadVec(DataFile, MeasPosYName, PosY);
        ReadVec(DataFile, MeasPosZName, PosZ);
        assert(PosX.size() == PosY.size());
        assert(PosX.size() == PosZ.size());
      }

    void SaveScalarGravityMeasurements(const std::string &filename,
        ThreeDGravityModel::tScalarMeasVec &Data,
        ThreeDGravityModel::tMeasPosVec &PosX,
        ThreeDGravityModel::tMeasPosVec &PosY,
        ThreeDGravityModel::tMeasPosVec &PosZ)
      {
        assert(Data.size() == PosX.size());
        assert(Data.size() == PosY.size());
        assert(Data.size() == PosZ.size());

        NcFile DataFile(filename.c_str(), NcFile::Replace);

        NcDim *StatNumDim;
        MakeStatNumDim(Data.size(),StatNumDim,DataFile);

        WriteVec(DataFile, MeasPosXName, PosX, StatNumDim);
        WriteVec(DataFile, MeasPosYName, PosY, StatNumDim);
        WriteVec(DataFile, MeasPosYName, PosY, StatNumDim);

        //Write the measurements
        NcVar *DataVar = DataFile.add_var(ScalarGravityName.c_str(), ncDouble,
            StatNumDim);
        DataVar->add_att("units", "m/s^2");
        DataVar->add_att("_FillValue", -1.0);

        DataVar->put(&Data[0], StatNumDim->size());
      }

    void ReadScalarGravityMeasurements(const std::string &filename,
        ThreeDGravityModel::tScalarMeasVec &Data,
        ThreeDGravityModel::tMeasPosVec &PosX,
        ThreeDGravityModel::tMeasPosVec &PosY,
        ThreeDGravityModel::tMeasPosVec &PosZ)
      {
        NcFile DataFile(filename.c_str(), NcFile::ReadOnly);
        ReadVec(DataFile, MeasPosXName, PosX);
        ReadVec(DataFile, MeasPosYName, PosY);
        ReadVec(DataFile, MeasPosZName, PosZ);
        ReadVec(DataFile, ScalarGravityName, Data);
      }

    void ReadTensorGravityMeasurements(const std::string &filename,
        ThreeDGravityModel::tTensorMeasVec &Data,
        ThreeDGravityModel::tMeasPosVec &PosX,
        ThreeDGravityModel::tMeasPosVec &PosY,
        ThreeDGravityModel::tMeasPosVec &PosZ)
      {
        NcFile DataFile(filename.c_str(), NcFile::ReadOnly);
        ReadVec(DataFile, MeasPosXName, PosX);
        ReadVec(DataFile, MeasPosYName, PosY);
        ReadVec(DataFile, MeasPosZName, PosZ);
        assert(PosX.size() == PosY.size());
        assert(PosX.size() == PosZ.size());

        rmat GravMat(3, 3);
        Data.assign(PosX.size(), GravMat);
        ReadMatComp(DataFile, UxxName, Data, 0, 0);
        ReadMatComp(DataFile, UxyName, Data, 0, 1);
        ReadMatComp(DataFile, UxzName, Data, 0, 2);
        ReadMatComp(DataFile, UyxName, Data, 1, 0);
        ReadMatComp(DataFile, UyyName, Data, 1, 1);
        ReadMatComp(DataFile, UyzName, Data, 1, 2);
        ReadMatComp(DataFile, UzxName, Data, 2, 0);
        ReadMatComp(DataFile, UzyName, Data, 2, 1);
        ReadMatComp(DataFile, UzzName, Data, 2, 2);
      }

    void SaveTensorGravityMeasurements(const std::string &filename,
        ThreeDGravityModel::tTensorMeasVec &Data,
        ThreeDGravityModel::tMeasPosVec &PosX,
        ThreeDGravityModel::tMeasPosVec &PosY,
        ThreeDGravityModel::tMeasPosVec &PosZ)
      {
        assert(PosX.size() == PosY.size());
        assert(PosX.size() == PosZ.size());
        assert(PosX.size() == Data.size());

        NcFile DataFile(filename.c_str(), NcFile::Replace);
        NcDim *StatNumDim;
        MakeStatNumDim(Data.size(),StatNumDim,DataFile);
        WriteVec(DataFile, MeasPosXName, PosX, StatNumDim);
        WriteVec(DataFile, MeasPosYName, PosY, StatNumDim);
        WriteVec(DataFile, MeasPosZName, PosZ, StatNumDim);

        WriteMatComp(DataFile, UxxName, Data, 0, 0, StatNumDim);
        WriteMatComp(DataFile, UxyName, Data, 0, 1, StatNumDim);
        WriteMatComp(DataFile, UxzName, Data, 0, 2, StatNumDim);
        WriteMatComp(DataFile, UyxName, Data, 1, 0, StatNumDim);
        WriteMatComp(DataFile, UyyName, Data, 1, 1, StatNumDim);
        WriteMatComp(DataFile, UyzName, Data, 1, 2, StatNumDim);
        WriteMatComp(DataFile, UzxName, Data, 2, 0, StatNumDim);
        WriteMatComp(DataFile, UzyName, Data, 2, 1, StatNumDim);
        WriteMatComp(DataFile, UzzName, Data, 2, 2, StatNumDim);
      }
  }
