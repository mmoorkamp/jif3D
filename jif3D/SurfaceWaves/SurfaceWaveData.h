/*
 * SurfaceWaveData.h
 *
 *  Created on: 19 Sep 2019
 *      Author: bweise
 */

#ifndef SURFACEWAVES_SURFACEWAVEDATA_H_
#define SURFACEWAVES_SURFACEWAVEDATA_H_

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
      std::vector<double> GetStatComb() const
        {
          return stat_comb;
        }
      std::vector<double> GetEventStatComb() const
        {
          return event_stat_comb;
        }
      double GetDummy() const
        {
          return dummy;
        }
      double GetCentrLon() const
        {
          return lon_centr;
        }
      virtual void WriteNetCDF(const std::string &filename) override;
      void WriteStationLocations(const std::string &filename);
    private:
      // Earthquake locations
      std::vector<double> EventPosX, EventPosY, EventPosZ;
      // Periods
      std::vector<double> periods;
      // Event-station-combinations
      std::vector<double> stat_comb, event_stat_comb;
      double lon_centr, dummy;
      };
  }


#endif /* SURFACEWAVES_SURFACEWAVEDATA_H_ */
