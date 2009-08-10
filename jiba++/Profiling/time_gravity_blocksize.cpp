#include <omp.h>
#include <iostream>
#include <fstream>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/assign/std/vector.hpp>
#include "../Global/convert.h"
#include "../Gravity/ThreeDGravityModel.h"
#include "../Gravity/ThreeDGravityCalculator.h"
#include "../Gravity/FullSensitivityGravityCalculator.h"
#include "../Gravity/MinMemGravityCalculator.h"
#include "../Gravity/ThreeDGravityFactory.h"

using namespace boost::assign;

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
            Model.SetDensities()[i][j][k] = 0.1 + double(rand() % 1000) / 300.0;
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


    const size_t nrunspersize = 5;
    std::string filename = "blocksize.time";
    boost::shared_ptr<jiba::ThreeDGravityCalculator> Calculator;
    boost::shared_ptr<jiba::TensorCudaGravityImp> Implementation(
        new jiba::TensorCudaGravityImp);

    std::vector<double> blocksizes;

    blocksizes += 64, 76, 90, 107, 128, 152, 181, 192, 215, 255;
    const size_t nruns = blocksizes.size();
    std::ofstream outfile(filename.c_str());
    std::cout << " Starting calculations. " << std::endl;
    for (size_t i = 0; i < nruns; ++i)
      {
        const size_t modelsize = 60;
        std::cout << "Blocksize: " << blocksizes.at(i) << std::endl;
        jiba::ThreeDGravityModel GravityTest;
        Implementation->SetCUDABlockSize(blocksizes.at(i));
        Calculator = boost::shared_ptr<jiba::ThreeDGravityCalculator>(
            new jiba::MinMemGravityCalculator(Implementation));
        double rawruntime = 0.0;
        for (size_t j = 0; j < nrunspersize; ++j)
          {
            rawruntime = 0.0;
            MakeTestModel(GravityTest, modelsize);

            boost::posix_time::ptime firststarttime =
                boost::posix_time::microsec_clock::local_time();
            jiba::rvec gravmeas(Calculator->Calculate(GravityTest));

            boost::posix_time::ptime firstendtime =
                boost::posix_time::microsec_clock::local_time();
            rawruntime += (firstendtime - firststarttime).total_microseconds();
          }
        rawruntime /= nrunspersize;
        outfile << blocksizes.at(i) << " " << rawruntime << std::endl;
      }

  }
