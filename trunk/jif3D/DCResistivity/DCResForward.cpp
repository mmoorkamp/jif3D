#include <iostream>
#include <string>
#include "../Global/FileUtil.h"
#include "../Global/Noise.h"
#include "ThreeDDCResistivityModel.h"
#include "DCResistivityCalculator.h"
#include "ReadWriteDCResistivityData.h"

int main()
  {
    std::string ModelFilename = jif3D::AskFilename("DCResistivity Model File: ", true);

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

    jif3D::DCResistivityCalculator Calculator;
    jif3D::rvec SynthData = Calculator.Calculate(Model);
    double error = 0.0;
    std::cout << "DC Apparent Resistivity error (s): ";
    std::cin >> error;
    //if we want to add noise to the data
    jif3D::rvec Errors(SynthData.size(), 0.0);
    if (error > 0.0)
      {
        jif3D::AddNoise(SynthData, 0.0, error);
        std::fill(Errors.begin(), Errors.end(), error);
      }
    //std::string DataFilename = jif3D::AskFilename("Output file: ", false);
    jif3D::SaveApparentResistivity(ModelFilename +".dcdata.nc", SynthData, Errors, Model);
    Model.WriteVTK(ModelFilename + ".vtk");

  }
