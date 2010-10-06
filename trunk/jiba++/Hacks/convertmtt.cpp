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
        jiba::rvec Impedances, Errors;
        std::vector<double> Frequencies, StatXCoord(nstats),
            StatYCoord(nstats), StatZCoord(nstats);
        StationFile.open(filename.c_str());
        std::vector<double> CurrFrequencies;
        jiba::rvec CurrImpedances, CurrErrors;
        std::string StationName;
        StationFile >> StatXCoord.front() >> StatYCoord.front()
            >> StatZCoord.front() >> StationName;
        if (StationFile.good())
          {
            jiba::ReadImpedancesFromMTT(StationName, Frequencies,
                CurrImpedances, CurrErrors);
          }

        const size_t nfreq = Frequencies.size();
        assert(nfreq * 8 == CurrImpedances.size());
        Impedances.resize(nstats * nfreq * 8);
        Errors.resize(nstats * nfreq * 8);
        size_t stationindex = 0;
        std::cout << stationindex << " " << StationName << " "
            << CurrFrequencies.size() << " " << nfreq << std::endl;
        while (StationFile.good() && stationindex < nstats)
          {
            for (size_t i = 0; i < nfreq; ++i)
              {
                std::copy(CurrImpedances.begin() + i * 8,
                    CurrImpedances.begin() + (i + 1) * 8, Impedances.begin()
                        + i * nstats * 8 + stationindex * 8);
                std::copy(CurrErrors.begin() + i * 8, CurrErrors.begin() + (i
                    + 1) * 8, Errors.begin() + i * nstats * 8 + stationindex
                    * 8);
              }
            double xcoord, ycoord, zcoord;
            StationFile >> xcoord >> ycoord >> zcoord >> StationName;

            if (StationFile.good())
              {
                StatXCoord.at(stationindex + 1) = xcoord;
                StatYCoord.at(stationindex + 1) = ycoord;
                StatZCoord.at(stationindex + 1) = zcoord;
                jiba::ReadImpedancesFromMTT(StationName, CurrFrequencies,
                    CurrImpedances, CurrErrors);
                std::cout << stationindex + 1 << " " << StationName << " "
                    << CurrFrequencies.size() << " " << nfreq << std::endl;
                if (CurrFrequencies.size() != nfreq)
                  throw jiba::FatalException(
                      "Number of frequencies in current file does not match number of frequencies in first file !");

              }
            ++stationindex;
          }
        std::string outfilename = jiba::AskFilename("Output file: ", false);
        jiba::WriteImpedancesToNetCDF(outfilename, Frequencies, StatXCoord,
            StatYCoord, StatZCoord, Impedances, Errors);
      }
  }
