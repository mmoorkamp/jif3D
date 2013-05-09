//============================================================================
// Name        : tomoregion.cpp
// Author      : 27 Mar 2012
// Version     : 
// Copyright   : 2012, mm489
//============================================================================

#include "../Tomo/ReadWriteTomographyData.h"
#include "../Tomo/ThreeDSeismicModel.h"
#include "../Global/FileUtil.h"
#include <iostream>

int main()
  {
    std::string datafilename = jiba::AskFilename("Data Filename: ");

    jiba::rvec TomoData;
    jiba::ThreeDSeismicModel Orig, New;
    jiba::ReadTraveltimes(datafilename, TomoData, Orig);

    double minx, maxx, miny, maxy;
    std::cout << "Minimum Northing value: ";
    std::cin >> minx;
    std::cout << "Maximum Northing value: ";
    std::cin >> maxx;
    std::cout << "Minimum Easting value: ";
    std::cin >> miny;
    std::cout << "Maximum Easting value: ";
    std::cin >> maxy;

    std::vector<size_t> NewMeasIndex(Orig.GetMeasPosX().size()), NewSourceIndex(
        Orig.GetSourcePosX().size());
    size_t currindex = 0;
    for (size_t i = 0; i < Orig.GetMeasPosX().size(); ++i)
      {
        if (Orig.GetMeasPosX().at(i) < maxx && Orig.GetMeasPosX().at(i) > minx
            && Orig.GetMeasPosY().at(i) < maxy && Orig.GetMeasPosY().at(i) > miny)
          {
            New.AddMeasurementPoint(Orig.GetMeasPosX().at(i), Orig.GetMeasPosY().at(i),
                Orig.GetMeasPosZ().at(i));
            NewMeasIndex.at(i) = currindex;
            ++currindex;
          }
      }

    currindex = 0;
    for (size_t i = 0; i < Orig.GetSourcePosX().size(); ++i)
      {
        if (Orig.GetSourcePosX().at(i) < maxx && Orig.GetSourcePosX().at(i) > minx
            && Orig.GetSourcePosY().at(i) < maxy && Orig.GetSourcePosY().at(i) > miny)
          {
            New.AddSource(Orig.GetSourcePosX().at(i), Orig.GetSourcePosY().at(i),
                Orig.GetSourcePosZ().at(i));
            NewSourceIndex.at(i) = currindex;
            ++currindex;
          }
      }

    size_t nmeas = TomoData.size();
    std::vector<double> TmpData;
    for (size_t i = 0; i < nmeas; ++i)
      {
        const size_t OldRecIndex = Orig.GetReceiverIndices().at(i);
        const size_t OldSourceIndex = Orig.GetSourceIndices().at(i);
        if (Orig.GetMeasPosX().at(OldRecIndex) < maxx
            && Orig.GetMeasPosX().at(OldRecIndex) > minx
            && Orig.GetMeasPosY().at(OldRecIndex) < maxy
            && Orig.GetMeasPosY().at(OldRecIndex) > miny
            && Orig.GetSourcePosX().at(OldSourceIndex) < maxx
            && Orig.GetSourcePosX().at(OldSourceIndex) > minx
            && Orig.GetSourcePosY().at(OldSourceIndex) < maxy
            && Orig.GetSourcePosY().at(OldSourceIndex) > miny)
          {
            TmpData.push_back(TomoData(i));
            New.AddMeasurementConfiguration(NewSourceIndex.at(OldSourceIndex),NewMeasIndex.at(OldRecIndex));
          }
      }
    jiba::rvec NewTomoData(TmpData.size());
    std::copy(TmpData.begin(),TmpData.end(),NewTomoData.begin());
    std::string newdatafilename = datafilename + ".cut.nc";
    jiba::SaveTraveltimes(newdatafilename,NewTomoData,New);

  }

