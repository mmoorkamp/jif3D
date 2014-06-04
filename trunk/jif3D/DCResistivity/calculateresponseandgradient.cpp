//============================================================================
// Name        : calculateresponseandgradient.cpp
// Author      : May, 2014
// Version     :
// Copyright   : 2014, zhanjie shi
//============================================================================
#include <iostream>
#include <string>
#include <cassert>
#include <fstream>
#include "../Global/FileUtil.h"
#include "ThreeDDCResistivityModel.h"
#include "DCResistivityCalculator.h"
#include "ReadWriteDCResistivityData.h"

int main()
  {
    std::string ModelFilename = jif3D::AskFilename("Model NetCDF File: ", true);
    jif3D::ThreeDDCResistivityModel Model;
    Model.ReadNetCDF(ModelFilename);

    std::string insourcefilename = jif3D::AskFilename("Source ASCII file: ", false);
    std::ifstream insourcefile(insourcefilename.c_str());
    std::string inreceiverfilename = jif3D::AskFilename("Receiver ASCII file: ", false);
    std::ifstream inreceiverfile(inreceiverfilename.c_str());
    double sposposx, sposposy, sposposz, snegposx, snegposy, snegposz;
    double receiverposx1, receiverposy1, receiverposz1, receiverposx2, receiverposy2,
        receiverposz2;
    int sourceindics;
    while (insourcefile.good())
      {
        insourcefile >> sposposx >> sposposy >> sposposz >> snegposx >> snegposy
            >> snegposz;
        if (insourcefile.good())
          {
            Model.AddSource(sposposx, sposposy, sposposz, snegposx, snegposy, snegposz);
          }
      }
    assert(Model.GetSourcePosPosX().size() == Model.GetSourcePosPosY().size());
    assert(Model.GetSourcePosPosX().size() == Model.GetSourcePosPosZ().size());
    assert(Model.GetSourceNegPosX().size() == Model.GetSourceNegPosY().size());
    assert(Model.GetSourceNegPosX().size() == Model.GetSourceNegPosZ().size());
    assert(Model.GetSourcePosPosX().size() == Model.GetSourceNegPosX().size());
    while (inreceiverfile.good())
      {
        inreceiverfile >> receiverposx1 >> receiverposy1 >> receiverposz1 >> receiverposx2
            >> receiverposy2 >> receiverposz2 >> sourceindics;
        if (inreceiverfile.good())
          {
            Model.AddMeasurementPoint(receiverposx1, receiverposy1, receiverposz1,
                receiverposx2, receiverposy2, receiverposz2, sourceindics);
          }
      }
    assert(Model.GetMeasPosX().size() == Model.GetMeasPosY().size());
    assert(Model.GetMeasPosX().size() == Model.GetMeasPosZ().size());
    assert(Model.GetMeasSecPosX().size() == Model.GetMeasSecPosY().size());
    assert(Model.GetMeasSecPosX().size() == Model.GetMeasSecPosZ().size());
    assert(Model.GetMeasPosX().size() == Model.GetMeasSecPosX().size());
    assert(Model.GetMeasPosX().size() == Model.GetSourceIndices().size());

    std::string DataFilename = jif3D::AskFilename("Response NetCDF file: ", false);
    std::string outresponseASCIIfilename = jif3D::AskFilename("Response ASCII file: ",
        false);
    std::string MisfitFilename = jif3D::AskFilename("Misfit ASCII file: ", false);
    std::string GradientASCIIfilename = jif3D::AskFilename("Gradient ASCII file: ",
        false);

    jif3D::DCResistivityCalculator Calculator;
    jif3D::rvec SynthData = Calculator.Calculate(Model);
    jif3D::rvec Error(SynthData.size(), 0.0);
    jif3D::SaveApparentResistivity(DataFilename, SynthData, Error, Model);
    std::ofstream outresponseASCIIfile(outresponseASCIIfilename.c_str());
    for (int i = 0; i < SynthData.size(); i++)
      {
        outresponseASCIIfile << SynthData[i] << "\n";
      }

    std::ifstream inmisfitfile(MisfitFilename.c_str());
    std::vector<double> misfitdata;
    double misfitdatatemp;
    while (inmisfitfile.good())
      {
        inmisfitfile >> misfitdatatemp;
        if (inmisfitfile.good())
          {
            misfitdata.push_back(misfitdatatemp);
          }
      }
    jif3D::rvec Misfit(misfitdata.size(), 0.0);
    for (int i = 0; i < misfitdata.size(); i++)
      {
        Misfit[i] = misfitdata[i];
      }
    jif3D::DCResistivityCalculator Gradientcalculator;
    jif3D::rvec GradientData = Gradientcalculator.LQDerivative(Model, Misfit);
    std::ofstream outgradientASCIIfile(GradientASCIIfilename.c_str());
    for (int i = 0; i < GradientData.size(); i++)
      {
        outgradientASCIIfile << GradientData[i] << "\n";
      }

  }
