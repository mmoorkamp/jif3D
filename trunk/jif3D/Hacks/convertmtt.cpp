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
#include "../MT/MTData.h"
#include "../MT/TipperData.h"
int main()
  {
    std::string filename = jif3D::AskFilename("File with station names/positions: ");
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
    std::string WantDist, DistFilename;
    std::cout << "Read in Distortion File: ";
    std::cin >> WantDist;
    std::vector<double> C;
    if (WantDist != "n")
      {
        DistFilename = jif3D::AskFilename("Distortion File: ");
        std::ifstream distfile(DistFilename.c_str());
        std::copy(std::istream_iterator<double>(distfile),
            std::istream_iterator<double>(), std::back_inserter(C));
        if (C.size() != nstats * 4)
          {
            std::cerr << "Number of elements in Distortion file " << C.size()
                << " does not match 4* number of stations " << nstats << std::endl;
            return 100;
          }
      }
    else
      {
        C.resize(nstats * 4, 0);
        for (size_t i = 0; i < nstats; ++i)
          {
            C.at(i * 4) = 1.0;
            C.at(i * 4 + 3) = 1.0;
          }
      }

    std::cout << "Number of stations: " << nstats << std::endl;
    if (nstats > 0)
      {
        std::vector<double> Impedances, Errors, Tipper, TipErr;
        std::vector<double> Frequencies, StatXCoord(nstats), StatYCoord(nstats),
            StatZCoord(nstats);
        std::vector<std::string> Names;
        StationFile.open(filename.c_str());
        std::vector<double> CurrFrequencies;
        std::vector<double> CurrImpedances, CurrErrors, CurrTip, CurrTipErr;
        std::string StationName;
        StationFile >> StatXCoord.front() >> StatYCoord.front() >> StatZCoord.front()
            >> StationName;
        if (StationFile.good())
          {
            std::string extension = jif3D::GetFileExtension(StationName);
            if (extension == ".mtt")
              {
                jif3D::ReadImpedancesFromMTT(StationName, Frequencies, CurrImpedances,
                    CurrErrors, CurrTip, CurrTipErr);
              }
            else
              {
                double XC, YC, ZC;
                jif3D::ReadImpedancesFromJ(StationName, Frequencies, XC, YC, ZC,
                    CurrImpedances, CurrErrors, CurrTip, CurrTipErr);
              }
          }
        Names.push_back(StationName);
        const size_t nfreq = Frequencies.size();
        assert(nfreq * 8 == CurrImpedances.size());
        Impedances.resize(nstats * nfreq * 8);
        Errors.resize(nstats * nfreq * 8);
        Tipper.resize(nstats * nfreq * 4, 0.0);
        TipErr.resize(nstats * nfreq * 4, 1.0);
        size_t stationindex = 0;
        std::cout << stationindex << " " << StationName << " " << CurrFrequencies.size()
            << " " << nfreq << std::endl;
        while (StationFile.good() && stationindex < nstats)
          {
            for (size_t i = 0; i < nfreq; ++i)
              {
                std::copy(CurrImpedances.begin() + i * 8,
                    CurrImpedances.begin() + (i + 1) * 8,
                    Impedances.begin() + i * nstats * 8 + stationindex * 8);
                std::copy(CurrErrors.begin() + i * 8, CurrErrors.begin() + (i + 1) * 8,
                    Errors.begin() + i * nstats * 8 + stationindex * 8);

                if (!CurrTip.empty())
                  {
                    std::copy(CurrTip.begin() + i * 4, CurrTip.begin() + (i + 1) * 4,
                        Tipper.begin() + i * nstats * 4 + stationindex * 4);
                    std::copy(CurrTipErr.begin() + i * 4,
                        CurrTipErr.begin() + (i + 1) * 4,
                        TipErr.begin() + i * nstats * 4 + stationindex * 4);
                  }
              }
            double xcoord, ycoord, zcoord;
            StationFile >> xcoord >> ycoord >> zcoord >> StationName;

            if (StationFile.good())
              {
                StatXCoord.at(stationindex + 1) = xcoord;
                StatYCoord.at(stationindex + 1) = ycoord;
                StatZCoord.at(stationindex + 1) = zcoord;
                Names.push_back(StationName);
                std::string extension = jif3D::GetFileExtension(StationName);
                if (extension == ".mtt")
                  {
                    jif3D::ReadImpedancesFromMTT(StationName, CurrFrequencies,
                        CurrImpedances, CurrErrors, CurrTip, CurrTipErr);
                  }
                else
                  {
                    double XC, YC, ZC;
                    jif3D::ReadImpedancesFromJ(StationName, CurrFrequencies, XC, YC, ZC,
                        CurrImpedances, CurrErrors, CurrTip, CurrTipErr);
                  }
                std::cout << stationindex + 1 << " " << StationName << " "
                    << CurrFrequencies.size() << " " << nfreq << std::endl;
                if (CurrFrequencies.size() != nfreq)
                  throw jif3D::FatalException(
                      "Number of frequencies in current file does not match number of frequencies in first file !");

              }
            ++stationindex;
          }
        std::string outfilename = jif3D::AskFilename("Output file: ", false);
        jif3D::MTData MTData;
        MTData.SetMeasurementPoints(StatXCoord, StatYCoord, StatZCoord);
        MTData.SetDataAndErrors(Impedances, Errors);
        MTData.SetFrequencies(Frequencies);
        MTData.SetNames(Names);
        MTData.CompleteObject();
        MTData.WriteNetCDF(outfilename);
        MTData.WriteMeasurementPoints(outfilename+".vtk");

        jif3D::TipperData TipperData;
        TipperData.SetMeasurementPoints(StatXCoord, StatYCoord, StatZCoord);
        TipperData.SetDataAndErrors(Tipper, TipErr);
        TipperData.SetFrequencies(Frequencies);
        TipperData.CompleteObject();
        TipperData.WriteNetCDF(outfilename + ".tip.nc");

      }
  }
