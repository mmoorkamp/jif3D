#include <omp.h>
#include <iostream>
#include <fstream>
#include "../Global/convert.h"
#include "../Gravity/ThreeDGravityModel.h"
#include "../Gravity/ThreeDGravityCalculator.h"
#include "../Gravity/FullSensitivityGravityCalculator.h"
#include "../Gravity/MinMemGravityCalculator.h"
#include "../Gravity/ThreeDGravityFactory.h"

void MakeTestModel(jiba::ThreeDGravityModel &Model, const size_t size)
  {
    Model.SetXCellSizes().resize(boost::extents[size]);
    Model.SetYCellSizes().resize(boost::extents[size]);
    Model.SetZCellSizes().resize(boost::extents[size]);

    for (size_t i = 0; i < size; ++i) // set the values of the inner cells
      {
        Model.SetXCellSizes()[i] = rand() % 10000 + 1000;
        Model.SetYCellSizes()[i] = rand() % 10000 + 1000;
        Model.SetZCellSizes()[i] = 500;
      }
    Model.SetDensities().resize(boost::extents[size][size][size]);
    for (size_t i = 0; i < size; ++i)
      for (size_t j = 0; j < size; ++j)
        for (size_t k = 0; k < size; ++k)
          {
            Model.SetDensities()[i][j][k] = 1.0 + double(rand() % 1000) / 300.0;
          }
    const size_t nmeas = 30;
    for (size_t i = 0; i < nmeas; ++i)
      Model.AddMeasurementPoint(rand() % 50000 + 2e4, rand() % 50000 + 2e4, 0.0);
    std::vector<double> bg_dens, bg_thick;
    bg_dens.push_back(1.0);
    bg_dens.push_back(1.0);
    bg_dens.push_back(5.0);
    bg_dens.push_back(5.0);
    bg_thick.push_back(200.0);
    bg_thick.push_back(300.0);
    bg_thick.push_back(3500.0);
    bg_thick.push_back(1000.0);
    Model.SetBackgroundDensities(bg_dens);
    Model.SetBackgroundThicknesses(bg_thick);
  }

int main(int ac, char* av[])
  {

    const size_t nruns = 5;

    std::string filename = ("ftgprec.out");

    boost::shared_ptr<jiba::ThreeDGravityCalculator>
        GPUCalculator = jiba::CreateGravityCalculator<
            jiba::MinMemGravityCalculator>::MakeTensor(true);

  //  boost::shared_ptr<jiba::ThreeDGravityCalculator>
   //     CPUCalculator = jiba::CreateGravityCalculator<
    //        jiba::MinMemGravityCalculator>::MakeTensor(false);

    std::ofstream outfile(filename.c_str());
    std::cout << " Starting calculations. " << std::endl;
    for (size_t i = 0; i < nruns; ++i)
      {
        const size_t modelsize = 80; //30 + rand() % 50;
        std::cout << "Current model size: " << pow(modelsize, 3) << std::endl;
        jiba::ThreeDGravityModel GravityTest;

        MakeTestModel(GravityTest, modelsize);

        jiba::rvec gpugravmeas(GPUCalculator->Calculate(GravityTest));
       // jiba::rvec cpugravmeas(CPUCalculator->Calculate(GravityTest));
        //for (size_t j = 0; j < gpugravmeas.size(); ++j)
         // {
          //  outfile << (gpugravmeas(j) - cpugravmeas(j)) / cpugravmeas(j)
           //     << std::endl;
         // }
      }

  }

