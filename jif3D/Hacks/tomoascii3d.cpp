//============================================================================
// Name        : picksort.cpp
// Author      : Sep 28, 2010
// Version     : 
// Copyright   : 2010, mmoorkamp
// Modified    : July 2017, mpaulatto
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
    std::vector<double> TTError;

//    while (PickFile.good())
//      {
//        double CurrShotIndex, CurrRecIndex, CurrTT;
//        PickFile >> CurrShotIndex >> CurrRecIndex;
//        char Line[255];
//        PickFile.getline(Line, 255);
//        try
//          {
//            jif3D::convert(Line, CurrTT);
//            if (PickFile.good())
//              {
//                if (CurrTT > 0.0)
//                  {
//                    ShotIndex.push_back(floor(CurrShotIndex));
//                    RecIndex.push_back(floor(CurrRecIndex));
//                    TravelTime.push_back(CurrTT);
//                  }
//             }
//          } catch (...)
//          {
//            //here we catch the occurence of NaN in the file
//          }
//      }

    while (PickFile.good())
      {
        double CurrShotIndex, CurrRecIndex, CurrTT, CurrDt;
        PickFile >> CurrShotIndex >> CurrRecIndex >> CurrTT >> CurrDt;
        if (PickFile.good())
          {
             ShotIndex.push_back(floor(CurrShotIndex));
             RecIndex.push_back(floor(CurrRecIndex));
             TravelTime.push_back(CurrTT);
             TTError.push_back(CurrDt);
          }
       }



    typedef std::map<size_t, std::vector<double> > MyMap;
    MyMap SourceMap, RecMap;

    const size_t nshot = ShotNo.size();
    for (size_t i = 0; i < nshot; ++i)
      {
        std::vector<double> SourcePos(3, 0.0);
        SourcePos[0] = SourceX.at(i);
        SourcePos[1] = SourceY.at(i);
        SourcePos[2] = SourceZ.at(i);
        SourceMap.insert(std::make_pair(ShotNo.at(i), SourcePos));
      }

    const size_t nrec = RecNo.size();
    for (size_t i = 0; i < nrec; ++i)
      {
        std::vector<double> RecPos(3, 0.0);
        RecPos[0] = RecX.at(i);
        RecPos[1] = RecY.at(i);
        RecPos[2] = RecZ.at(i);
        RecMap.insert(std::make_pair(RecNo.at(i), RecPos));
      }

    jif3D::ThreeDSeismicModel Model;
//    const double depth = 12.0;

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
        double offset = std::sqrt(
            std::pow(SourceX.at(SI) - RecX.at(RI), 2)
                + std::pow(SourceY.at(SI) - RecY.at(RI), 2));
        double distance = std::sqrt(
            std::pow(SourceX.at(SI) - RecX.at(RI), 2)
                + std::pow(SourceY.at(SI) - RecY.at(RI), 2)
                + std::pow(SourceZ.at(SI) - RecZ.at(RI), 2));
        double currvel = distance / TravelTime.at(i);
        if (offset > mindist && currvel > minvel && currvel < maxvel)
          {
            Model.AddMeasurementPoint(RecX.at(RI), RecY.at(RI), RecZ.at(RI));
            Model.AddMeasurementConfiguration(SI, measindex);
            ++measindex;
          }
        else
          {
            TravelTime.at(i) = -1.0;
            TTError.at(i) = - 1.0;
          }

      }
    TravelTime.erase(std::remove(TravelTime.begin(), TravelTime.end(), -1.0),
        TravelTime.end());
    TTError.erase(std::remove(TTError.begin(), TTError.end(), -1.0),
        TTError.end());
    std::cout << "NTimes in file: " << TravelTime.size() << std::endl;
// MP debug
//    std::cout << "Test checkpoint 1";
    std::string outfilename = jif3D::AskFilename("Output file: ", false);
// MP debug
//    std::cout << "Test checkpoint 2";
    jif3D::rvec TT(TravelTime.size());
    std::copy(TravelTime.begin(), TravelTime.end(), TT.begin());
//    jif3D::rvec Error(TT.size(), 0.0);
// This doesn't work if some of the traveltimes are culled because of the
// minimum offset or velocity. Inconsistency in array size.
    jif3D::rvec Error(TTError.size());
    std::copy(TTError.begin(), TTError.end(), Error.begin());
// MP debug
//    std::cout << Error.size();
//    std::cout << TT.size();
    jif3D::SaveTraveltimes(outfilename, TT, Error, Model);
  }


