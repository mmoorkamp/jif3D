/*
 * SurfaceWaveData.h
 *
 *  Created on: 19 Sep 2019
 *      Author: bweise
 */

#ifndef SURFACEWAVES_SURFACEWAVEDATA_H_
#define SURFACEWAVES_SURFACEWAVEDATA_H_

#include <map>
#include <tuple>
#include "../DataBase/GeneralData.h"
namespace jif3D
  {
    class SurfaceWaveData: public GeneralData
      {
    public:
      SurfaceWaveData();
      virtual void ReadNetCDF(const std::string &datafile) override;
      std::vector<double> GetEventPosLat() const
        {
          return EventPosLat;
        }
      std::vector<double> GetEventPosLon() const
        {
          return EventPosLon;
        }
      std::vector<double> GetEventPosZ() const
        {
          return EventPosZ;
        }
      std::vector<double> GetPeriods() const
        {
          return periods;
        }
      std::vector<int> GetStatPairs() const
        {
          return StationPairs;
        }
      std::multimap<int, std::tuple<int, int>> GetIndexMap() const
        {
          return indexmap;
        }
      std::vector<int> GetDataPerT() const
        {
          return NDataPerT;
        }
      double GetCentrLon() const
        {
          return lon_centr;
        }
      void SetEventPositions(const std::vector<double> &elat,
          const std::vector<double> &elon, const std::vector<double> &epz)
        {
          EventPosLat.resize(elat.size());
          EventPosLon.resize(elon.size());
          EventPosZ.resize(epz.size());
          std::copy(elat.begin(), elat.end(), EventPosLat.begin());
          std::copy(elon.begin(), elon.end(), EventPosLon.begin());
          std::copy(epz.begin(), epz.end(), EventPosZ.begin());
        }
      void SetPeriods(const std::vector<double> &T)
        {
          periods.resize(T.size());
          std::copy(T.begin(), T.end(), periods.begin());
        }
      void SetLonCentr(const double &lc)
        {
          lon_centr = lc;
        }
      void SetStationPairs(const std::vector<int> &station1,
          const std::vector<int> &station2)
        {
          StationPairs.resize(station1.size() * 2);
          std::copy(station1.begin(), station1.end(), StationPairs.begin());
          std::copy(station2.begin(), station2.end(),
              StationPairs.begin() + station1.size());
        }
      void SetNDataPerT(const std::vector<int> &ndata)
        {
          NDataPerT.resize(ndata.size());
          std::copy(ndata.begin(), ndata.end(), NDataPerT.begin());
        }
      void SetIndexMap(const std::vector<int> &PairInd, const std::vector<int> &EventInd,
          const std::vector<int> &PeriodInd)
        {
          indexmap.clear();
          for (size_t ndata = 0; ndata < PairInd.size(); ndata++)
            {
              indexmap.insert(
                  std::pair<int, std::tuple<int, int>>(PairInd[ndata],
                      std::make_tuple(EventInd[ndata], PeriodInd[ndata])));
            }
        }
      virtual void WriteNetCDF(const std::string &filename) override;
      void WriteStationLocations(const std::string &filename);
      void WriteRayPaths(const std::vector<std::vector<double>> &pathmap_e,
          const std::vector<std::vector<double>> &pathmap_n) const;
    private:
      // Earthquake locations
      std::vector<double> EventPosLat, EventPosLon, EventPosZ;
      // Periods
      std::vector<double> periods;
      // Event-station-combinations
      std::vector<int> StationPairs, NDataPerT;
      std::multimap<int, std::tuple<int, int>> indexmap;
      double lon_centr;
      };
  }

#endif /* SURFACEWAVES_SURFACEWAVEDATA_H_ */
