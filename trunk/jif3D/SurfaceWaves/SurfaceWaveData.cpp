/*
 * SurfaceWaveData.cpp
 *
 *  Created on: 19 Sep 2019
 *      Author: bweise
 */

#include "../SurfaceWaves/SurfaceWaveData.h"
#include "../ModelBase/VTKTools.h"
#include "../Global/NetCDFPortHelper.h"
#include "../Global/NetCDFTools.h"
#include <netcdf>
#include <vector>
#include <tuple>
#include <map>
using netCDF::NcFile;
using netCDF::NcVar;
using netCDF::NcVarAtt;
using netCDF::NcDim;

namespace jif3D
  {

    void SurfaceWaveData::ReadNetCDF(const std::string &datafile)
      {

        // Read phase delay time observations
        NcFile dtpFile(datafile, NcFile::read);

        ReadVec(dtpFile, "Periods", periods);

        std::vector<double> MPX, MPY, MPZ;
        ReadVec(dtpFile, "MeasPosX", MPX);
        ReadVec(dtpFile, "MeasPosY", MPY);
        ReadVec(dtpFile, "MeasPosZ", MPZ);
        SetMeasurementPoints(MPX, MPY, MPZ);

        ReadVec(dtpFile, "EventPosX", EventPosX);
        ReadVec(dtpFile, "EventPosY", EventPosY);
        ReadVec(dtpFile, "EventPosY", EventPosZ);
        ReadVec(dtpFile, "TraveltimesPerPeriod", NDataPerT);

        NcDim nsrcsIn = dtpFile.getDim("NumberOfPairs");
        NcDim SRIn = dtpFile.getDim("SR");
        StationPairs.resize(nsrcsIn.getSize() * SRIn.getSize());
        NcVar src_rcvr_cmbIn = dtpFile.getVar("StationPairs");
        src_rcvr_cmbIn.getVar(StationPairs.data());

        std::vector<int> PairInd, TInd, EventInd;
        ReadVec(dtpFile, "PairIndex", PairInd);
        ReadVec(dtpFile, "PeriodIndex", TInd);
        ReadVec(dtpFile, "EventIndex", EventInd);

        std::vector<double> dtp;
        ReadVec(dtpFile, "dtp", dtp);

        std::vector<double> err;
        std::string errname = "dtp_error";
        err.resize(dtp.size(), 0.0);
        // check if there is an error in the data file, set to zero if not.
        try
          {
            if (!dtpFile.getVar(errname.c_str()).isNull())
              {
                NcVar errIn = dtpFile.getVar("dtau");
                errIn.getVar(err.data());
              }
          } catch (netCDF::exceptions::NcException &ex)
          {
            // ignore
          }

        // Organize data in multimap
        for (int ndata = 0; ndata < dtp.size(); ndata++)
          {
            datamap.insert(
                std::pair<int, std::tuple<int, int, double, double>>(PairInd[ndata],
                    std::make_tuple(EventInd[ndata], TInd[ndata], dtp[ndata],
                        err[ndata])));
          }

        std::vector<double> dtp_sorted(dtp.size()), err_sorted(dtp.size());
        std::multimap<int, std::tuple<int, int, double, double>>::iterator it;
        for (it = datamap.begin(); it != datamap.end(); ++it)
          {
            auto datatuple = (*it).second;
            double tmpdat = std::get<2>(datatuple);
            double tmperr = std::get<3>(datatuple);
            dtp_sorted.push_back(tmpdat);
            err_sorted.push_back(tmperr);
          }

        SetDataAndErrors(dtp_sorted, err_sorted);

        NcVar mpnIn = dtpFile.getVar("MeasPosX");
        NcVarAtt lonc = mpnIn.getAtt("Central_meridian");
        lonc.getValues(&lon_centr);
      }

    SurfaceWaveData::SurfaceWaveData() :
        EventPosX(), EventPosY(), EventPosZ(), periods(), lon_centr(), NDataPerT(), StationPairs()
      {
      }

    void SurfaceWaveData::WriteNetCDF(const std::string &filename)
      {
        const size_t sr = 2;
        const size_t nevents = GetEventPosX().size();
        const size_t nperiods = GetPeriods().size();
        const size_t nstats = GetMeasPosX().size();
        const size_t npairs = GetStatPairs().size();
        const size_t ntimes = GetData().size();
        const std::string lon_centr = "Central_meridian";

        //create a netcdf file
        NcFile DataFile(filename, NcFile::replace);
        //we use the station number as a dimension
        NcDim NumberOfEvents = DataFile.addDim("NumberOfEvents", nevents);
        NcDim NumberOfPairs = DataFile.addDim("NumberOfPairs", npairs);
        NcDim NumberOfTraveltimes = DataFile.addDim("NumberOfTraveltimes", ntimes);
        NcDim NumberOfStations = DataFile.addDim("NumberOfStations", nstats);
        NcDim NumberOfPeriods = DataFile.addDim("NumberOfPeriods", nperiods);
        NcDim SR = DataFile.addDim("SR", sr);

        //write out the event coordinates
        WriteVec(DataFile, "EventPosX", GetEventPosX(), NumberOfEvents, "m");
        WriteVec(DataFile, "EventPosY", GetEventPosY(), NumberOfEvents, "m");
        WriteVec(DataFile, "EventPosZ", GetEventPosZ(), NumberOfEvents, "m");
        //write out the positions of the stations, i.e. measurement positions
        NcVar MeasPosXVar = DataFile.addVar("MeasPosX", netCDF::ncDouble,
            NumberOfStations);
        MeasPosXVar.putAtt("units", "m");
        MeasPosXVar.putAtt(lon_centr, NC_DOUBLE, GetCentrLon());
        jif3D::cxxport::put_legacy_ncvar(MeasPosXVar, GetMeasPosX().data(), nstats);

        NcVar MeasPosYVar = DataFile.addVar("MeasPosY", netCDF::ncDouble,
            NumberOfStations);
        MeasPosYVar.putAtt("units", "m");
        MeasPosYVar.putAtt(lon_centr, NC_DOUBLE, GetCentrLon());
        jif3D::cxxport::put_legacy_ncvar(MeasPosYVar, GetMeasPosY().data(), nstats);

        WriteVec(DataFile, "MeasPosZ", GetMeasPosZ(), NumberOfStations, "m");
        // write out periods
        WriteVec(DataFile, "Periods", GetPeriods(), NumberOfPeriods, "s");
        WriteVec(DataFile, "StationPairs", GetStatPairs(), NumberOfPairs, "");
        WriteVec(DataFile, "TraveltimesPerPeriod", GetDataPerT(), NumberOfPeriods, "");
        WriteVec(DataFile, "dtp", GetData(), NumberOfTraveltimes, "s");
        WriteVec(DataFile, "dtp_error", GetErrors(), NumberOfTraveltimes, "s");

        datamap = GetDataMap();
        std::vector<int> EInd(GetData().size()), PInd(GetData().size()), TInd(
            GetData().size());
        std::multimap<int, std::tuple<int, int, double, double>>::iterator it;
        for (it = datamap.begin(); it != datamap.end(); ++it)
          {
            int ptmp = (*it).first;
            auto datatuple = (*it).second;
            int etmp = std::get<0>(datatuple);
            int ttmp = std::get<1>(datatuple);
            PInd.push_back(ptmp);
            EInd.push_back(etmp);
            TInd.push_back(ttmp);
          }

        WriteVec(DataFile, "EventIndex", EInd, NumberOfTraveltimes, "");
        WriteVec(DataFile, "PairIndex", PInd, NumberOfTraveltimes, "");
        WriteVec(DataFile, "PeriodIndex", TInd, NumberOfTraveltimes, "");

      }

    void SurfaceWaveData::WriteStationLocations(const std::string &filename)
      {
        std::vector<double> SourceNum(GetMeasPosX().size());
        std::iota(SourceNum.begin(), SourceNum.end(), 1);
        jif3D::Write3DDataToVTK(filename, "Stations", SourceNum, GetMeasPosX(),
            GetMeasPosY(), GetMeasPosZ());
      }
  }
