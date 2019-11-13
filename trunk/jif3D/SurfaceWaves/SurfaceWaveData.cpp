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

        ReadVec(dtpFile,"Periods",periods);
        std::vector<double> MPX, MPY, MPZ;
        ReadVec(dtpFile,"MeasPosX",MPX);
        ReadVec(dtpFile,"MeasPosY",MPY);
        ReadVec(dtpFile,"MeasPosZ",MPZ);
        SetMeasurementPoints(MPX, MPY, MPZ);

        NcDim nsrcsIn = dtpFile.getDim("NumberOfRays");
        NcDim SRIn = dtpFile.getDim("SR");
        stat_comb.resize(nsrcsIn.getSize() * SRIn.getSize());
        NcVar src_rcvr_cmbIn=dtpFile.getVar("StatComb");
        src_rcvr_cmbIn.getVar(stat_comb.data());

        NcDim nperiodsIn = dtpFile.getDim("NumberOfPeriods");
        NcDim nevents_per_srcIn = dtpFile.getDim("EventsPerSRC");
        std::vector<double> dtp(nperiodsIn.getSize() * nevents_per_srcIn.getSize() * nsrcsIn.getSize()), err;
        NcVar dtpIn = dtpFile.getVar("dtp");
        dtpIn.getVar(dtp.data());

        std::string errname = "dtau";
        err.resize(dtp.size(),0.0);
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
        SetDataAndErrors(dtp, err);

        event_stat_comb.resize(nevents_per_srcIn.getSize() * nsrcsIn.getSize());
        NcVar event_stat_cmbIn=dtpFile.getVar("EventStatComb");
        event_stat_cmbIn.getVar(event_stat_comb.data());

        ReadVec(dtpFile,"EventPosX",EventPosX);
        ReadVec(dtpFile,"EventPosY",EventPosY);

        NcVarAtt dummyIn = dtpIn.getAtt("_FillValue");
        dummyIn.getValues(&dummy);
        
        NcVar mpnIn = dtpFile.getVar("MeasPosX");        
        NcVarAtt lonc = mpnIn.getAtt("Central_meridian");
        lonc.getValues(&lon_centr);
      }

    SurfaceWaveData::SurfaceWaveData()
      {
      }

    void SurfaceWaveData::WriteNetCDF(const std::string &filename)
      {
        const size_t sr = 2;
        const size_t nevents = GetEventPosX().size();
        const size_t nperiods = GetPeriods().size();
        const size_t nrays = GetStatComb().size() / sr;
        const size_t nstats = GetMeasPosX().size();
        const size_t events_per_src = GetEventStatComb().size() / nrays;
        const std::string dummy = "_FillValue";
        const std::string lon_centr = "Central_meridian";

        //create a netcdf file
        NcFile DataFile(filename, NcFile::replace);
        //we use the station number as a dimension
        NcDim NumberOfEvents = DataFile.addDim("NumberOfEvents", nevents);
        NcDim EventsPerSRC = DataFile.addDim("EventsPerSRC", events_per_src);
        NcDim NumberOfRays = DataFile.addDim("NumberOfRays", nrays);
        NcDim NumberOfStations = DataFile.addDim("NumberOfStations", nstats);
        NcDim NumberOfPeriods = DataFile.addDim("NumberOfPeriods", nperiods);
        NcDim SR = DataFile.addDim("SR", sr);

        //write out the event coordinates
        WriteVec(DataFile, "EventPosX", GetEventPosX(), NumberOfEvents, "m");
        WriteVec(DataFile, "EventPosY", GetEventPosY(), NumberOfEvents, "m");
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

        std::vector<NcDim> dimVec_esc;
        dimVec_esc.push_back(NumberOfRays);
        dimVec_esc.push_back(EventsPerSRC);
        NcVar EvStaCoVar = DataFile.addVar("EventStatComb", netCDF::ncDouble,
            dimVec_esc);
        EvStaCoVar.putAtt(dummy, NC_DOUBLE, GetDummy());
        cxxport::put_legacy_ncvar(EvStaCoVar, GetEventStatComb().data(), nrays, events_per_src);

        std::vector<NcDim> dimVec_sc;
        dimVec_sc.push_back(NumberOfRays);
        dimVec_sc.push_back(SR);
        NcVar StaCoVar = DataFile.addVar("StatComb", netCDF::ncDouble, dimVec_sc);
        cxxport::put_legacy_ncvar(StaCoVar, GetStatComb().data(), nrays, sr);

        std::vector<NcDim> dimVec_dtp;
        dimVec_dtp.push_back(NumberOfRays);
        dimVec_dtp.push_back(EventsPerSRC);
        dimVec_dtp.push_back(NumberOfPeriods);
        NcVar DTPVar = DataFile.addVar("dtp", netCDF::ncDouble, dimVec_dtp);
        DTPVar.putAtt(dummy, NC_DOUBLE, GetDummy());
        cxxport::put_legacy_ncvar(DTPVar, GetData().data(), nrays, events_per_src, nperiods);
      }

    void SurfaceWaveData::WriteStationLocations(const std::string &filename)
      {
        std::vector<double> SourceNum(GetMeasPosX().size());
        std::iota(SourceNum.begin(), SourceNum.end(), 1);
        jif3D::Write3DDataToVTK(filename, "Stations", SourceNum, GetMeasPosX(),
            GetMeasPosY(), GetMeasPosZ());
      }
  }
