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
    std::vector<size_t> ShotNo, RecNo;
    std::vector<double> SourceX, SourceY, SourceZ, RecX, RecY, RecZ;

    std::string RecFileName = jif3D::AskFilename("File with receiver positions: ");
    std::ifstream RecFile(RecFileName.c_str());
    while (RecFile.good())
      {
        size_t CurrRecNo;
        double CurrRecX, CurrRecY, CurrRecZ;
        RecFile >> CurrRecNo >> CurrRecX >> CurrRecY >> CurrRecZ;
        if (RecFile.good())
          {
            RecNo.push_back(CurrRecNo);
            RecX.push_back(CurrRecX);
            RecY.push_back(CurrRecY);
            RecZ.push_back(CurrRecZ);
          }
      }

    std::string SorFileName = jif3D::AskFilename("File with source positions: ");
    std::ifstream SorFile(SorFileName.c_str());
    while (SorFile.good())
      {
        size_t CurrSorNo;
        double CurrSorX, CurrSorY, CurrSorZ;
        SorFile >> CurrSorNo >> CurrSorX >> CurrSorY >> CurrSorZ;
        if (SorFile.good())
          {
            ShotNo.push_back(CurrSorNo);
            SourceX.push_back(CurrSorX);
            SourceY.push_back(CurrSorY);
            SourceZ.push_back(CurrSorZ);
          }
      }

    std::string PickFileName = jif3D::AskFilename("Pick-file name: ");
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
            jif3D::convert(Line, CurrTT);
            if (PickFile.good())
              {
                if (CurrTT > 0.0)
                  {
                    ShotIndex.push_back(floor(CurrShotIndex));
                    RecIndex.push_back(floor(CurrRecIndex));
                    TravelTime.push_back(CurrTT);
                  }

              }
          } catch (...)
          {
            //here we catch the occurence of NaN in the file
          }
      }

    typedef std::map<size_t, std::vector<double> > MyMap;
    MyMap SourceMap, RecMap;

    const size_t nshot = ShotNo.size();
    for (size_t i = 0; i < nshot; ++i)
      {
        std::vector<double> SourcePos({SourceX.at(i), SourceY.at(i), SourceZ.at(i)});
        SourceMap.insert(std::make_pair(ShotNo.at(i), SourcePos));
      }

    const size_t nrec = RecNo.size();
    for (size_t i = 0; i < nrec; ++i)
      {
        std::vector<double> RecPos({RecX.at(i),RecY.at(i),RecZ.at(i)});
        RecMap.insert(std::make_pair(RecNo.at(i), RecPos));
      }

    jif3D::ThreeDSeismicModel Model;
    for (auto source : SourceMap)
      {
        Model.AddSource(source.second[0], source.second[1], source.second[2]);
      }



    const size_t ntime = TravelTime.size();
    std::cout << "NTimes total: " << ntime << std::endl;
    size_t measindex = 0;
    double mindist = 0.0;
    double minvel = 0.0;
    const double maxvel = 10000;
    std::cout << "Minimum offset [m]: ";
    std::cin >> mindist;
    std::cout << "Minimum velocity [m/s]: ";
    std::cin >> minvel;
    for (size_t i = 0; i < ntime; ++i)
      {
        MyMap::iterator siter = SourceMap.find(ShotIndex.at(i));
        MyMap::iterator riter = RecMap.find(RecIndex.at(i));
        size_t SI = std::distance(SourceMap.begin(),siter);
        size_t RI = std::distance(RecMap.begin(),riter);
        double distance = std::sqrt(
            std::pow(SourceX.at(SI) - RecX.at(RI), 2)
                + std::pow(SourceY.at(SI) - RecY.at(RI), 2));
        double currvel = distance / TravelTime.at(i);
        if (distance > mindist && currvel > minvel && currvel < maxvel)
          {
            Model.AddMeasurementPoint(RecX.at(RI), RecY.at(RI), RecZ.at(RI));
            Model.AddMeasurementConfiguration(SI, measindex);
            ++measindex;
          }
        else
          {
            TravelTime.at(i) = -1.0;
          }

      }
    TravelTime.erase(std::remove(TravelTime.begin(), TravelTime.end(), -1.0),
        TravelTime.end());
    std::cout << "NTimes in file: " << TravelTime.size() << std::endl;
    std::string outfilename = jif3D::AskFilename("Output file: ", false);
    jif3D::rvec TT(TravelTime.size());
    std::copy(TravelTime.begin(), TravelTime.end(), TT.begin());
    jif3D::rvec Error(TT.size(), 0.0);
    jif3D::SaveTraveltimes(outfilename, TT, Error, Model);
  }

