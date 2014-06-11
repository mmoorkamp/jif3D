#include <iostream>
#include <string>
#include "../Global/FileUtil.h"
#include "ThreeDDCResistivityModel.h"
#include "DCResistivityCalculator.h"
#include "ReadWriteDCResistivityData.h"

int main()
  {
    std::string ModelFilename = jif3D::AskFilename("Model File: ", true);

    jif3D::ThreeDDCResistivityModel Model;
    Model.ReadNetCDF(ModelFilename);

    size_t nelecx, nelecy;

    std::cout << "N Electrodes in x-direction: ";
    std::cin >> nelecx;
    std::cout << "N Electrodes in y-direction: ";
    std::cin >> nelecy;

    double minx, miny, deltax, deltay;
    std::cout << "First position in x-direction: ";
    std::cin >> minx;
    std::cout << "First position in y-direction: ";
    std::cin >> miny;

    std::cout << "Spacing in x-direction: ";
    std::cin >> deltax;
    std::cout << "Spacing in y-direction: ";
    std::cin >> deltay;

    size_t nsources;
    std::cout << "N Sources: ";
    std::cin >> nsources;

    std::vector<int> SourceIndices;
    for (size_t i = 0; i < nsources; ++i)
      {
        double index;
        std::cout << "Index " << i << " : ";
        std::cin >> index;
        SourceIndices.push_back(index);
      }

    for (size_t i = 0; i < nsources; ++i)
      {
        size_t Index = SourceIndices[i];
        double xpos = minx + deltax * (Index % nelecx);
        double ypos = miny + deltay * (Index / nelecy);
        Model.AddSource(xpos, ypos, 0.0, xpos + deltax, ypos, 0.0);
        for (size_t j = 0; j < nelecx * nelecy; ++j)
          {
            if (j != Index)
              {
                double measposx = minx + deltax * (j % nelecx);
                double measposy = miny + deltay * (j / nelecy);
                Model.AddMeasurementPoint(measposx, measposy, 0.0, measposx + deltax,
                    measposy + deltay, 0.0, Index);
              }
          }
      }

    jif3D::DCResistivityCalculator Calculator;

    jif3D::rvec SynthData = Calculator.Calculate(Model);
    jif3D::rvec Error(SynthData.size(), 0.0);

    std::string DataFilename = jif3D::AskFilename("Output file: ", false);
    jif3D::SaveApparentResistivity(DataFilename, SynthData, Error, Model);

  }
