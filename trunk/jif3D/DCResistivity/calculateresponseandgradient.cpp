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
#include "../DCResistivity/DCResistivityData.h"

int main()
  {
    std::string ModelFilename = jif3D::AskFilename("Model NetCDF File: ", true);
    jif3D::ThreeDDCResistivityModel DCModel;
    jif3D::DCResistivityData DCData;
    DCModel.ReadNetCDF(ModelFilename);

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
        	DCData.AddSource(sposposx, sposposy, sposposz, snegposx, snegposy, snegposz);
          }
      }
    assert(DCData.GetSourcePosPosX().size() == DCData.GetSourcePosPosY().size());
    assert(DCData.GetSourcePosPosX().size() == DCData.GetSourcePosPosZ().size());
    assert(DCData.GetSourceNegPosX().size() == DCData.GetSourceNegPosY().size());
    assert(DCData.GetSourceNegPosX().size() == DCData.GetSourceNegPosZ().size());
    assert(DCData.GetSourcePosPosX().size() == DCData.GetSourceNegPosX().size());
    while (inreceiverfile.good())
      {
        inreceiverfile >> receiverposx1 >> receiverposy1 >> receiverposz1 >> receiverposx2
            >> receiverposy2 >> receiverposz2 >> sourceindics;
        if (inreceiverfile.good())
          {
        	DCData.AddMeasurementPoint(receiverposx1, receiverposy1, receiverposz1,
                receiverposx2, receiverposy2, receiverposz2, sourceindics);
          }
      }
    assert(DCData.GetMeasPosX().size() == DCData.GetMeasPosY().size());
    assert(DCData.GetMeasPosX().size() == DCData.GetMeasPosZ().size());
    assert(DCData.GetMeasSecPosX().size() == DCData.GetMeasSecPosY().size());
    assert(DCData.GetMeasSecPosX().size() == DCData.GetMeasSecPosZ().size());
    assert(DCData.GetMeasPosX().size() == DCData.GetMeasSecPosX().size());
    assert(DCData.GetMeasPosX().size() == DCData.GetSourceIndices().size());

    std::string DataFilename = jif3D::AskFilename("Response NetCDF file: ", false);
    std::string outresponseASCIIfilename = jif3D::AskFilename("Response ASCII file: ",
        false);
    std::string MisfitFilename = jif3D::AskFilename("Misfit ASCII file: ", false);
    std::string GradientASCIIfilename = jif3D::AskFilename("Gradient ASCII file: ",
        false);

    jif3D::DCResistivityCalculator Calculator;
    jif3D::rvec SynthData = Calculator.Calculate(DCModel, DCData);
    jif3D::rvec Error(SynthData.size(), 0.0);

    DCData.SetDataAndErrors(std::vector<double>(SynthData.begin(), SynthData.end()),
        std::vector<double>(Error.begin(), Error.end()));
    DCData.WriteNetCDF(DataFilename + ".dcdata.nc");
    std::ofstream outresponseASCIIfile(outresponseASCIIfilename.c_str());
    for (double data : SynthData)
      {
        outresponseASCIIfile << data << "\n";
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
    std::copy(misfitdata.begin(), misfitdata.end(), Misfit.begin());

    jif3D::DCResistivityCalculator Gradientcalculator;
    jif3D::rvec GradientData = Gradientcalculator.LQDerivative(DCModel, DCData, Misfit);
    std::ofstream outgradientASCIIfile(GradientASCIIfilename.c_str());
    for (double gradient : GradientData)
      {
        outgradientASCIIfile << gradient << "\n";
      }

  }
