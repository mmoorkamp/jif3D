//============================================================================
// Name        : ReadWriteMagneticData.cpp
// Author      : Max Moorkamp
// Version     :
// Copyright   : 2008, MM
//============================================================================

#include <netcdf>
#include "ReadWriteMagneticData.h"
#include "../Global/NumUtil.h"
#include "../Global/NetCDFTools.h"
#include "../ModelBase/NetCDFModelTools.h"

using netCDF::NcFile;
using netCDF::NcVar;
using netCDF::NcDim;

namespace jif3D
  {

    //! Read a component of the magnetic field for all measurement positions
    void ReadMagComp(NcFile &NetCDFFile, const std::string &CompName, rvec &MatVec,
        size_t n)
      {
        //read in the component from the netcdf file
        rvec tempdata(MatVec.size() / 3);
        ReadVec(NetCDFFile, CompName, tempdata);
        //copy to the right location in the matrix
        for (size_t i = 0; i < tempdata.size(); ++i)
          MatVec(i * 3 + n) = tempdata(i);
      }

    //! Write one component of the magnetic field for all measurement positions
    void WriteMagComp(NcFile &NetCDFFile, const std::string &CompName, const rvec &MatVec,
        const size_t n, NcDim &Dimension)
      {
        rvec tempdata(MatVec.size() / 3);
        for (size_t i = 0; i < tempdata.size(); ++i)
          {
            tempdata(i) = MatVec(i * 3 + n);
          }

        WriteVec(NetCDFFile, CompName, tempdata, Dimension, "1/s2");
      }

    static const std::string BxName = "Bx";
    static const std::string BxErrorName = "dBx";
    static const std::string ByName = "By";
    static const std::string ByErrorName = "dBy";
    static const std::string BzName = "Bz";
    static const std::string BzErrorName = "dBz";
    static const std::string TotalFieldName = "T";
    static const std::string TotalFieldErrorName = "dT";
    static const std::vector<std::string> ComponentNames =
      { BxName, ByName, BzName };
    static const std::vector<std::string> ComponentErrorNames =
      { BxErrorName, ByErrorName, BzErrorName };

    /*! Write the scalar Magnetic measurements and their position to a netcdf file
     * @param filename The name of the file including ending
     * @param Data The vector of scalar Magnetic measurements in m/s2
     * @param PosX The x-coordinate (Northing) of each measurement in m
     * @param PosY The y-coordinate (Easting) of each measurement in m
     * @param PosZ The z-coordinate (Depth) of each measurement in m
     * @param Error The error estimate for each measurement
     */
    void SaveTotalFieldMagneticMeasurements(const std::string &filename,
        const std::vector<double> &Data, const std::vector<double> &PosX,
        const std::vector<double> &PosY, const std::vector<double> &PosZ,
        const std::vector<double> &Error)
      {
        //make sure all vectors have consistent sizes
        assert(Data.size() == PosX.size());
        assert(Data.size() == PosY.size());
        assert(Data.size() == PosZ.size());

        //create a netcdf file
        NcFile DataFile(filename, NcFile::replace);
        //we use the station number as a dimension
        NcDim StatNumDim = DataFile.addDim(StationNumberName, Data.size());

        //write out the measurement coordinates
        WriteVec(DataFile, MeasPosXName, PosX, StatNumDim, "m");
        WriteVec(DataFile, MeasPosYName, PosY, StatNumDim, "m");
        WriteVec(DataFile, MeasPosZName, PosZ, StatNumDim, "m");
        //Write the measurements and the error
        WriteVec(DataFile, TotalFieldName, Data, StatNumDim, "nT");
        WriteVec(DataFile, TotalFieldErrorName, Error, StatNumDim, "nT");
      }

    /*! Read scalar Magnetic measurements and their position from a netcdf file
     * @param filename The name of the file including ending
     * @param Data The vector of scalar Magnetic measurements in m/s2
     * @param PosX The x-coordinate (Northing) of each measurement in m
     * @param PosY The y-coordinate (Easting) of each measurement in m
     * @param PosZ The z-coordinate (Depth) of each measurement in m
     * @param Error The error estimate for each measurement
     */
    void ReadTotalFieldMagneticMeasurements(const std::string &filename,
        std::vector<double> &Data, std::vector<double> &PosX, std::vector<double> &PosY,
        std::vector<double> &PosZ, std::vector<double> &Error)
      {
        NcFile DataFile(filename, NcFile::read);

        ReadVec(DataFile, MeasPosXName, PosX);
        ReadVec(DataFile, MeasPosYName, PosY);
        ReadVec(DataFile, MeasPosZName, PosZ);
        ReadVec(DataFile, TotalFieldName, Data);

        try
          {
            if (!DataFile.getVar(TotalFieldErrorName).isNull())
              {
                ReadVec(DataFile, TotalFieldErrorName, Error);
              }
            else
              {
                Error.resize(Data.size());
                Error.clear();
              }
          } catch (netCDF::exceptions::NcException &ex)
          {
            // ignore
          }
      }

    /*! Read each magnetic field component and their position to a netcdf file. Data will have
     * 3 consecutive entries (the three field components Bx, By, Bz ) per element in PosX,PosY and PosZ.
     * @param filename The name of the file including ending
     * @param Data The vector of magnetic measurements in nT
     * @param PosX The x-coordinate (Northing) of each measurement in m
     * @param PosY The y-coordinate (Easting) of each measurement in m
     * @param PosZ The z-coordinate (Depth) of each measurement in m
     * @param Error The error estimate for each measurement
     */
    void ReadMagneticComponentMeasurements(const std::string &filename, jif3D::rvec &Data,
        std::vector<double> &PosX, std::vector<double> &PosY, std::vector<double> &PosZ,
        jif3D::rvec &Error)
      {
        NcFile DataFile(filename, NcFile::read);

        ReadVec(DataFile, MeasPosXName, PosX);
        ReadVec(DataFile, MeasPosYName, PosY);
        ReadVec(DataFile, MeasPosZName, PosZ);

        assert(PosX.size() == PosY.size());
        assert(PosX.size() == PosZ.size());

        Data.resize(PosX.size() * 3);
        Error.resize(Data.size());
        Error.clear();
        for (size_t i = 0; i < 3; ++i)
          {
            ReadMagComp(DataFile, ComponentNames.at(i), Data, i);
          }
        try
          {
            for (size_t i = 0; i < 3; ++i)
              {
                if (!DataFile.getVar(ComponentErrorNames.at(i)).isNull())
                  {
                    ReadMagComp(DataFile, ComponentErrorNames.at(i), Error, i);
                  }
              }
          } catch (const netCDF::exceptions::NcException &ex)
          {
            // ignore
          }
      }

    /*! Write each magnetic field component and their position to a netcdf file. Data has to have
     * 3 consecutive entries (the three field components) per element in PosX,PosY and PosZ.
     * @param filename The name of the file including ending
     * @param Data The vector of magnetic measurements in nt
     * @param PosX The x-coordinate (Northing) of each measurement in m
     * @param PosY The y-coordinate (Easting) of each measurement in m
     * @param PosZ The z-coordinate (Depth) of each measurement in m
     * @param Error The error estimate for each measurement
     */
    void SaveTensorMagneticMeasurements(const std::string &filename,
        const jif3D::rvec &Data, const std::vector<double> &PosX,
        const std::vector<double> &PosY, const std::vector<double> &PosZ,
        const jif3D::rvec &Error)
      {
        const size_t nmeas = PosX.size();

        assert(nmeas == PosY.size());
        assert(nmeas == PosZ.size());
        assert(nmeas * 3 == Data.size());

        NcFile DataFile(filename, NcFile::replace);
        //create a dimension for the stations
        NcDim StatNumDim = DataFile.addDim(StationNumberName.c_str(), nmeas);
        //write out the measurement positions
        WriteVec(DataFile, MeasPosXName, PosX, StatNumDim, "m");
        WriteVec(DataFile, MeasPosYName, PosY, StatNumDim, "m");
        WriteVec(DataFile, MeasPosZName, PosZ, StatNumDim, "m");

        //write the three components of the Magnetic field
        for (size_t i = 0; i < 3; ++i)
          {
            WriteMagComp(DataFile, ComponentNames.at(i), Data, i, StatNumDim);
          }
        for (size_t i = 0; i < 3; ++i)
          {
            WriteMagComp(DataFile, ComponentErrorNames.at(i), Error, i, StatNumDim);
          }

      }
  }
