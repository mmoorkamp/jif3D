//============================================================================
// Name        : picksort.cpp
// Author      : Sep 28, 2010
// Version     : 
// Copyright   : 2010, mmoorkamp
//============================================================================

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <set>
#include <boost/math/special_functions/round.hpp>

#include "../Global/VecMat.h"
#include "../Global/FileUtil.h"
#include "../Global/convert.h"
#include "../Tomo/ThreeDSeismicModel.h"
#include "../Tomo/ReadWriteTomographyData.h"

size_t MakeKey(size_t RecNo, size_t SourceNo)
  {
    return 1e6 * RecNo + SourceNo;
  }

void ExtractIndex(size_t key, size_t &RecNo, size_t &SourceNo)
  {
    RecNo = key / 1e6;
    SourceNo = key - 1e6 * RecNo;
  }

int main()
  {
    std::string PosFileName = jiba::AskFilename("File with positions: ");
    std::ifstream PosFile(PosFileName.c_str());
    std::vector<size_t> ShotNo, RecNo, Tracl;
    std::vector<double> SourceX, SourceY, RecX, RecY;
    while (PosFile.good())
      {
        size_t CurrShotNo, CurrRecNo, CurrTracl;
        double CurrSourceX, CurrSourceY, CurrRecX, CurrRecY;
        PosFile >> CurrShotNo >> CurrTracl >> CurrRecNo >> CurrSourceY
            >> CurrSourceX >> CurrRecY >> CurrRecX;
        if (PosFile.good())
          {
            ShotNo.push_back(CurrShotNo);
            RecNo.push_back(CurrRecNo);
            Tracl.push_back(CurrTracl);
            SourceX.push_back(CurrSourceX);
            SourceY.push_back(CurrSourceY);
            RecX.push_back(CurrRecX);
            RecY.push_back(CurrRecY);
          }
      }

    typedef std::map<size_t, std::vector<double> > MyMap;
    MyMap SourceRecMap, SourceMap;

    const size_t nrecords = ShotNo.size();
    for (size_t i = 0; i < nrecords; ++i)
      {
        std::vector<double> SourceRecPos(4, 0.0), SourcePos(2, 0.0);
        SourceRecPos[0] = RecX.at(i);
        SourceRecPos[1] = RecY.at(i);
        SourceRecPos[2] = SourceX.at(i);
        SourceRecPos[3] = SourceY.at(i);
        SourcePos[0] = SourceX.at(i);
        SourcePos[1] = SourceY.at(i);
        size_t key = MakeKey(RecNo.at(i), ShotNo.at(i));
        SourceRecMap.insert(std::make_pair(key, SourceRecPos));
        SourceMap.insert(std::make_pair(ShotNo.at(i), SourcePos));
      }

    std::string PickFileName = jiba::AskFilename("Pick-file name: ");
    std::ifstream PickFile(PickFileName.c_str());
    std::vector<size_t> ShotIndex, RecIndex;
    std::vector<double> TravelTime;
    while (PickFile.good())
      {
        double CurrShotIndex, CurrRecIndex, CurrTT;
        PickFile >> CurrShotIndex >> CurrRecIndex;
        char Line[255];
        PickFile.getline(Line, 255);
        try
          {
            jiba::convert(Line, CurrTT);
            if (PickFile.good())
              {

                ShotIndex.push_back(floor(CurrShotIndex));
                RecIndex.push_back(floor(CurrRecIndex));
                TravelTime.push_back(CurrTT);

              }
          } catch (...)
          {
            //here we catch the occurence of NaN in the file
          }
      }

    jiba::ThreeDSeismicModel Model;
    size_t nrec = SourceRecMap.size();
    const double depth = 10.0;

    for (MyMap::iterator sourceiter = SourceMap.begin(); sourceiter
        != SourceMap.end(); ++sourceiter)
      {
        Model.AddSource(sourceiter->second[0], sourceiter->second[1], depth);
      }

    const size_t ntime = TravelTime.size();
    std::cout << "NTimes total: " << ntime << std::endl;
    size_t measindex = 0;
    double mindist = 0.0;
    std::cout << "Minimum offset: ";
    std::cin >> mindist;
    for (size_t i = 0; i < ntime; ++i)
      {
        std::cout << RecIndex.at(i) << " " << ShotIndex.at(i) << std::endl;
        const size_t Key = MakeKey(RecIndex.at(i), ShotIndex.at(i));
        MyMap::iterator sriter = SourceRecMap.find(Key);

        if (sriter != SourceRecMap.end())
          {

            MyMap::iterator sourceiter = SourceMap.find(ShotIndex.at(i));

            size_t CurrShotIndex = std::distance(SourceMap.begin(), sourceiter);
            if (sourceiter->second[0] != sriter->second[2])
              {
                std::cerr << "Mismatch in X-component of source: "
                    << CurrShotIndex << " " << sourceiter->second[0] << " "
                    << sriter->second[2] << std::endl;
                return 100;
              }
            if (sourceiter->second[1] != sriter->second[3])
              {
                std::cerr << "Mismatch in Y-component of source: "
                    << CurrShotIndex << " " << sourceiter->second[1] << " "
                    << sriter->second[3] << std::endl;
                return 100;
              }
            double distance = std::sqrt(std::pow(sourceiter->second[0]
                - sriter->second[0], 2) + std::pow(sourceiter->second[1]
                - sriter->second[1], 2));
            if (distance > mindist)
              {
                Model.AddMeasurementPoint(sriter->second[0], sriter->second[1],
                    depth);
                Model.AddMeasurementConfiguration(CurrShotIndex, measindex);
                ++measindex;
              }
            else
              {
                TravelTime.at(i) = -1.0;
              }
          }
        else
          {
            std::cout << "Could not find : ";
            std::cout << RecIndex.at(i) << " " << ShotIndex.at(i) << std::endl;
            TravelTime.at(i) = -1.0;
          }
      }
    TravelTime.erase(std::remove(TravelTime.begin(), TravelTime.end(), -1.0),
        TravelTime.end());
    std::cout << "NTimes in file: " << TravelTime.size() << std::endl;
    std::string outfilename = jiba::AskFilename("Output file: ", false);
    jiba::rvec TT(TravelTime.size());
    std::copy(TravelTime.begin(), TravelTime.end(), TT.begin());
    jiba::SaveTraveltimes(outfilename, TT, Model);
  }

