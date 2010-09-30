//============================================================================
// Name        : convertmtt.cpp
// Author      : Sep 21, 2010
// Version     : 
// Copyright   : 2010, mmoorkamp
//============================================================================

#include <iostream>
#include <fstream>
#include <string>
#include <boost/program_options.hpp>
#include "../Global/FileUtil.h"
#include "../Global/FatalException.h"
#include "../MT/ReadWriteImpedances.h"

int main()
  {
    std::string filename = jiba::AskFilename(
        "File with station names/positions: ");
    std::ifstream StationFile;

    StationFile.open(filename.c_str());
    char dummy[2048];
    size_t nstats = 0;
    while (StationFile.good())
      {
        StationFile.getline(dummy, 2048);
        if (StationFile.good())
          {
            ++nstats;
          }
      }
    StationFile.close();
    std::cout << "Number of stations: " << nstats << std::endl;
    if (nstats > 0)
      {
        jiba::rvec Impedances;
        std::vector<double> Frequencies, StatXCoord(nstats),
            StatYCoord(nstats), StatZCoord(nstats);
        StationFile.open(filename.c_str());
        std::vector<double> CurrFrequencies;
        jiba::rvec CurrImpedances;
        std::string StationName;
        StationFile >> StatXCoord.front() >> StatYCoord.front()
            >> StatZCoord.front() >> StationName;
        if (StationFile.good())
          {
            jiba::ReadImpedancesFromMTT(StationName, Frequencies,
                CurrImpedances);
          }

        const size_t nfreq = Frequencies.size();
        assert(nfreq * 8 == CurrImpedances.size());
        Impedances.resize(nstats * nfreq * 8);

        for (size_t i = 0; i < nfreq; ++i)
          {
            std::copy(CurrImpedances.begin() + i * 8, CurrImpedances.begin()
                + (i + 1) * 8, Impedances.begin() + i * nstats * 8);
          }
        size_t stationindex = 1;
        while (StationFile.good() && stationindex < nstats)
          {
            StationFile >> StatXCoord.at(stationindex) >> StatYCoord.at(
                stationindex) >> StatZCoord.at(stationindex) >> StationName;
            if (StationFile.good())
              {
                jiba::ReadImpedancesFromMTT(StationName, CurrFrequencies,
                    CurrImpedances);
                std::cout << stationindex << " " << StationName << " "
                    << CurrFrequencies.size() << " " << nfreq << std::endl;
                if (CurrFrequencies.size() != nfreq)
                  throw jiba::FatalException(
                      "Number of frequencies in current file does not match number of frequencies in first file !");

                for (size_t i = 0; i < nfreq; ++i)
                  {
                    std::copy(CurrImpedances.begin() + i * 8,
                        CurrImpedances.begin() + (i + 1) * 8,
                        Impedances.begin() + i * nstats * 8 + stationindex * 8);
                  }
              }
            ++stationindex;
          }
        std::string outfilename = jiba::AskFilename("Output file: ", false);
        jiba::WriteImpedancesToNetCDF(outfilename, Frequencies, StatXCoord,
            StatYCoord, StatZCoord, Impedances);
      }
  }
