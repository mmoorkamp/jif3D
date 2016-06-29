//============================================================================
// Name        : ReadWriteGravityData.cpp
// Author      : Max Moorkamp
// Version     :
// Copyright   : 2008, MM
//============================================================================

#include <netcdf>
#include "ReadWriteGravityData.h"
#include "../Global/NumUtil.h"
#include "../Global/NetCDFTools.h"
#include "../ModelBase/NetCDFModelTools.h"

using netCDF::NcFile;
using netCDF::NcDim;

namespace jif3D
  {

    static const std::string ScalarGravityName = "Scalar_gravity";
    static const std::string ScalarErrorName = "dGz";
//    static const std::string TensorErrorName = "dU";
    static const std::string UxxName = "Uxx";
    static const std::string UxyName = "Uxy";
    static const std::string UxzName = "Uxz";
    static const std::string UyxName = "Uyx";
    static const std::string UyyName = "Uyy";
    static const std::string UyzName = "Uyz";
    static const std::string UzxName = "Uzx";
    static const std::string UzyName = "Uzy";
    static const std::string UzzName = "Uzz";
    static const std::vector<std::string> TensorNames =
      { UxxName, UxyName, UxzName, UyxName, UyyName, UyzName, UzxName, UzyName, UzzName };

    //! Read a component of an FTG matrix for all measurement positions
    void ReadMatComp(NcFile &NetCDFFile, const std::string &CompName, rvec &MatVec,
        size_t n)
      {
        //read in the component from the netcdf file
        rvec tempdata(MatVec.size() / 9);
        ReadVec(NetCDFFile, CompName, tempdata);
        //copy to the right location in the matrix
        for (size_t i = 0; i < tempdata.size(); ++i)
          MatVec(i * 9 + n) = tempdata(i);
      }

    //! Write one component of an FTG matrix for all measurement positions
    void WriteMatComp(NcFile &NetCDFFile, const std::string &CompName, const rvec &MatVec,
        const size_t n, const NcDim &Dimension)
      {
        rvec tempdata(MatVec.size() / 9);

        for (size_t i = 0; i < tempdata.size(); ++i) {
          tempdata(i) = MatVec(i * 9 + n);
        }

        WriteVec(NetCDFFile, CompName, tempdata, Dimension, "1/s2");
      }

    /*! This function inspects the contents of a netcdf file to determine which kind
     * of gravity data is stored in it. If it encounters a variable called "Scalar_gravity"
     * it returns scalar, if it encounters a variable called "Uxx" it returns ftg. If none
     * of the two variables is encountered it returns unknown and if the file contains both,
     * tit returns whichever one it finds first.
     * @param filename The name of the netcdf file
     * @return One of scalar, ftg or unknown
     */
    GravityDataType IdentifyGravityDatafileType(const std::string &filename)
      {
        //open the file
        NcFile DataFile(filename.c_str(), NcFile::read);

        if(!DataFile.getVar(ScalarGravityName).isNull())
        {
          return scalar;
        }
        else if(!DataFile.getVar(UxxName).isNull())
        {
          return ftg;
        }


        return unknown;;
      }

    /*! Write the scalar gravity measurements and their position to a netcdf file
     * @param filename The name of the file including ending
     * @param Data The vector of scalar gravity measurements in m/s2
     * @param PosX The x-coordinate (Northing) of each measurement in m
     * @param PosY The y-coordinate (Easting) of each measurement in m
     * @param PosZ The z-coordinate (Depth) of each measurement in m
     * @param Error The error estimate for each measurement
     */
    void SaveScalarGravityMeasurements(const std::string &filename,
        const jif3D::rvec &Data, const ThreeDGravityModel::tMeasPosVec &PosX,
        const ThreeDGravityModel::tMeasPosVec &PosY,
        const ThreeDGravityModel::tMeasPosVec &PosZ, const jif3D::rvec &Error)
      {
        //make sure all vectors have consistent sizes
        assert(Data.size() == PosX.size());
        assert(Data.size() == PosY.size());
        assert(Data.size() == PosZ.size());
        jif3D::rvec LocalError(Error);
        if (Error.empty())
          {
            LocalError.resize(Data.size(), false);
            std::fill(LocalError.begin(), LocalError.end(), 0.0);
          }
        //create a netcdf file
        NcFile DataFile(filename.c_str(), NcFile::replace);
        //we use the station number as a dimension
        NcDim StatNumDim = DataFile.addDim(StationNumberName.c_str(), Data.size());

        //write out the measurement coordinates
        WriteVec(DataFile, MeasPosXName, PosX, StatNumDim, "m");
        WriteVec(DataFile, MeasPosYName, PosY, StatNumDim, "m");
        WriteVec(DataFile, MeasPosZName, PosZ, StatNumDim, "m");
        //Write the measurements and the error
        WriteVec(DataFile, ScalarGravityName, Data, StatNumDim, "m/s2");
        WriteVec(DataFile, ScalarErrorName, LocalError, StatNumDim, "m/s2");
      }

    /*! Read scalar gravity measurements and their position from a netcdf file
     * @param filename The name of the file including ending
     * @param Data The vector of scalar gravity measurements in m/s2
     * @param PosX The x-coordinate (Northing) of each measurement in m
     * @param PosY The y-coordinate (Easting) of each measurement in m
     * @param PosZ The z-coordinate (Depth) of each measurement in m
     * @param Error The error estimate for each measurement
     */
    void ReadScalarGravityMeasurements(const std::string &filename, jif3D::rvec &Data,
        ThreeDGravityModel::tMeasPosVec &PosX, ThreeDGravityModel::tMeasPosVec &PosY,
        ThreeDGravityModel::tMeasPosVec &PosZ, jif3D::rvec &Error)
      {
        NcFile DataFile(filename.c_str(), NcFile::read);
        ReadVec(DataFile, MeasPosXName, PosX);
        ReadVec(DataFile, MeasPosYName, PosY);
        ReadVec(DataFile, MeasPosZName, PosZ);
        ReadVec(DataFile, ScalarGravityName, Data);

        try {
          if (!DataFile.getVar(ScalarErrorName).isNull())
            {
              ReadVec(DataFile, ScalarErrorName, Error);
            }
          else
            {
              Error.resize(Data.size());
              Error.clear();
            }
        } catch(const netCDF::exceptions::NcException &ex) {
          // ignore
        }
      }

    /*! Read FTG measurements and their position from a netcdf file. Data will have
     * 9 consecutive entries (the nine components of the tensor) per element in PosX,PosY and PosZ.
     * @param filename The name of the file including ending
     * @param Data The vector of FTG measurements in 1/s2
     * @param PosX The x-coordinate (Northing) of each measurement in m
     * @param PosY The y-coordinate (Easting) of each measurement in m
     * @param PosZ The z-coordinate (Depth) of each measurement in m
     * @param Error The error estimate for each measurement
     */
    void ReadTensorGravityMeasurements(const std::string &filename, jif3D::rvec &Data,
        ThreeDGravityModel::tMeasPosVec &PosX, ThreeDGravityModel::tMeasPosVec &PosY,
        ThreeDGravityModel::tMeasPosVec &PosZ, jif3D::rvec &Error)
      {
        NcFile DataFile(filename.c_str(), NcFile::read);
        ReadVec(DataFile, MeasPosXName, PosX);
        ReadVec(DataFile, MeasPosYName, PosY);
        ReadVec(DataFile, MeasPosZName, PosZ);
        assert(PosX.size() == PosY.size());
        assert(PosX.size() == PosZ.size());

        Data.resize(PosX.size() * 9);
        Error.resize(Data.size());
        Error.clear();
        for (size_t i = 0; i < 9; ++i)
          {
            ReadMatComp(DataFile, TensorNames.at(i), Data, i);
          }

        for (size_t i = 0; i < 9; ++i)
          {
            try {
              const std::string currname = "d" + TensorNames.at(i);
              if (!DataFile.getVar(currname).isNull())
                {
                  ReadMatComp(DataFile, currname, Error, i);
                }
            } catch(const netCDF::exceptions::NcException &ex) {
              // ignore
            }
          }

      }

    /*! Write FTG measurements and their position to a netcdf file. Data has to have
     * 9 consecutive entries (the nine components of the tensor) per element in PosX,PosY and PosZ.
     * @param filename The name of the file including ending
     * @param Data The vector of FTG measurements in 1/s2
     * @param PosX The x-coordinate (Northing) of each measurement in m
     * @param PosY The y-coordinate (Easting) of each measurement in m
     * @param PosZ The z-coordinate (Depth) of each measurement in m
     * @param Error The error estimate for each measurement
     */
    void SaveTensorGravityMeasurements(const std::string &filename,
        const jif3D::rvec &Data, const ThreeDGravityModel::tMeasPosVec &PosX,
        const ThreeDGravityModel::tMeasPosVec &PosY,
        const ThreeDGravityModel::tMeasPosVec &PosZ, const jif3D::rvec &Error)
      {
        const size_t nmeas = PosX.size();
        assert(nmeas == PosY.size());
        assert(nmeas == PosZ.size());
        assert(nmeas * 9 == Data.size());

        NcFile DataFile(filename.c_str(), NcFile::replace);

        NcDim StatNumDim = DataFile.addDim(StationNumberName.c_str(), nmeas);

        WriteVec(DataFile, MeasPosXName, PosX, StatNumDim, "m");
        WriteVec(DataFile, MeasPosYName, PosY, StatNumDim, "m");
        WriteVec(DataFile, MeasPosZName, PosZ, StatNumDim, "m");
        //write the nine component of the gravity tensor
        for (size_t i = 0; i < 9; ++i)
          {
            WriteMatComp(DataFile, TensorNames.at(i), Data, i, StatNumDim);
          }
        for (size_t i = 0; i < 9; ++i)
          {
            WriteMatComp(DataFile, "d" + TensorNames.at(i), Error, i, StatNumDim);
          }

      }
  }
