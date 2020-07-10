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
      std::vector<double> GetEventPosX() const
        {
          return EventPosX;
        }
      std::vector<double> GetEventPosY() const
        {
          return EventPosY;
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
      std::multimap<int, std::tuple<int, int, double, double>> GetDataMap() const
        {
          return datamap;
        }
      std::vector<int> GetDataPerT() const
        {
          return NDataPerT;
        }
      double GetCentrLon() const
        {
          return lon_centr;
        }
      void SetEventPositions(const std::vector<double> &epx,
          const std::vector<double> &epy, const std::vector<double> &epz)
        {
          EventPosX.resize(epx.size());
          EventPosY.resize(epy.size());
          EventPosZ.resize(epz.size());
          std::copy(epx.begin(), epx.end(), EventPosX.begin());
          std::copy(epy.begin(), epy.end(), EventPosY.begin());
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
      void SetDataMap(const std::vector<int> &PairInd, const std::vector<int> &EventInd,
          const std::vector<int> &PeriodInd, const std::vector<double> &dtp,
          const std::vector<double> &error)
        {
          for (int ndata = 0; ndata < PairInd.size(); ndata++)
            {
              datamap.insert(
                  std::pair<int, std::tuple<int, int, double, double>>(PairInd[ndata],
                      std::make_tuple(EventInd[ndata], PeriodInd[ndata], dtp[ndata],
                          error[ndata])));
            }
          std::vector<double> dtp_sorted, err_sorted;
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
        }
      virtual void WriteNetCDF(const std::string &filename) override;
      void WriteStationLocations(const std::string &filename);
    private:
      // Earthquake locations
      std::vector<double> EventPosX, EventPosY, EventPosZ;
      // Periods
      std::vector<double> periods;
      // Event-station-combinations
      std::vector<int> StationPairs, NDataPerT;
      std::multimap<int, std::tuple<int, int, double, double>> datamap;
      double lon_centr;
      };
  }

#endif /* SURFACEWAVES_SURFACEWAVEDATA_H_ */
