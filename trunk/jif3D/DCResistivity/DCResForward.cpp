#include <iostream>
#include <string>
#include "../Global/FileUtil.h"
#include "../Global/Noise.h"
#include "ThreeDDCResistivityModel.h"
#include "DCResistivityCalculator.h"
#include "../DCResistivity/DCResistivityData.h"

int main()
  {
    std::string ModelFilename = jif3D::AskFilename("DCResistivity Model File: ", true);

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

    jif3D::DCResistivityCalculator Calculator;
    jif3D::rvec SynthData = Calculator.Calculate(DCModel, DCData);
    double error = 0.0;
    std::cout << "DC Apparent Resistivity error (s): ";
    std::cin >> error;
    //if we want to add noise to the data
    jif3D::rvec Errors(SynthData.size(), 0.0);
    if (error > 0.0)
      {
        jif3D::AddNoise(SynthData, 0.02, error);
        std::fill(Errors.begin(), Errors.end(), error);
      }
    //std::string DataFilename = jif3D::AskFilename("Output file: ", false);
    DCData.SetDataAndErrors(std::vector<double>(SynthData.begin(), SynthData.end()),
        std::vector<double>(Errors.begin(), Errors.end()));
    DCData.WriteNetCDF(ModelFilename + ".dcdata.nc");


    DCModel.WriteVTK(ModelFilename + ".vtk");

  }
