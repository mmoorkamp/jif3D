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
    SourceNo = key - 1e6*RecNo;
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
        PosFile >> CurrShotNo >> CurrTracl >> CurrRecNo >> CurrSourceX
            >> CurrSourceY >> CurrRecX >> CurrRecY;
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
    MyMap SourceRecMap;

    const size_t nrecords = ShotNo.size();
    for (size_t i = 0; i < nrecords; ++i)
      {
        std::vector<double> SourceRecPos(4, 0.0);
        SourceRecPos[0] = RecX.at(i);
        SourceRecPos[1] = RecY.at(i);
        SourceRecPos[2] = SourceX.at(i);
        SourceRecPos[3] = SourceY.at(i);
        size_t key = MakeKey(RecNo.at(i), ShotNo.at(i));
        SourceRecMap.insert(std::make_pair(key, SourceRecPos));
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

                ShotIndex.push_back(CurrShotIndex);
                RecIndex.push_back(CurrRecIndex);
                TravelTime.push_back(CurrTT);

              }
          } catch (...)
          {

          }
      }

    jiba::ThreeDSeismicModel Model;
    size_t nrec = SourceRecMap.size();

    std::cout << "Nrec: " << nrec << std::endl;
    for (MyMap::iterator RecIter = SourceRecMap.begin(); RecIter
        != SourceRecMap.end(); ++RecIter)
      {
        Model.AddMeasurementPoint(RecIter->second[0], RecIter->second[1], 10.0);
      }

    std::set<size_t> SourceSet;
    for (MyMap::iterator SourceIter = SourceRecMap.begin(); SourceIter
        != SourceRecMap.end(); ++SourceIter)
      {
        size_t RecIndex, SourceIndex;
        ExtractIndex(SourceIter->first, RecIndex, SourceIndex);
        if (SourceSet.insert(SourceIndex).second)
          {
            std::cout << "Inserted: " << SourceIndex << std::endl;
            Model.AddSource(SourceIter->second[2], SourceIter->second[3], 10.0);
          }
      }

    const size_t ntime = TravelTime.size();
    std::cout << "NTimes: " << ntime << std::endl;
    for (size_t i = 0; i < ntime; ++i)
      {
        const size_t Key = MakeKey(RecIndex.at(i), ShotIndex.at(i));
        MyMap::iterator sriter = SourceRecMap.find(Key);
        std::set<size_t>::iterator SourceIter = SourceSet.find(
            ShotIndex.at(i));
        if (sriter != SourceRecMap.end() && SourceIter != SourceSet.end())
          {
            size_t CurrRecIndex = std::distance(SourceRecMap.begin(), sriter);
            size_t CurrShotIndex = std::distance(SourceSet.begin(), SourceIter);
            Model.AddMeasurementConfiguration(CurrShotIndex, CurrRecIndex);
          }
        else
          {
            TravelTime.at(i) = -1.0;
          }
      }
    TravelTime.erase(std::remove(TravelTime.begin(),TravelTime.end(),-1.0),TravelTime.end());
    std::string outfilename = jiba::AskFilename("Output file: ", false);
    jiba::rvec TT(TravelTime.size());
    std::copy(TravelTime.begin(), TravelTime.end(), TT.begin());
    jiba::SaveTraveltimes(outfilename, TT, Model);
  }

