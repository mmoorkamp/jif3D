#include "../Global/Jif3DGlobal.h"
#include "../Global/convert.h"
#include "../Gravity/ThreeDGravityModel.h"
#include "../GravMag/ThreeDGravMagCalculator.h"
#include "../GravMag/FullSensitivityGravMagCalculator.h"
#include "../GravMag/MinMemGravMagCalculator.h"
#include "../Gravity/ThreeDGravityFactory.h"
#include <iostream>
#include <fstream>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/assign/std/vector.hpp>

using namespace boost::assign;



/*! \file time_gravity_blocksize.cpp
 * This program shows the execution time for the gavity forward calculation
 * using CUDA for a fixed model size, but different CUDA thread block sizes.
 */

void MakeTestModel(jif3D::ThreeDGravityModel &Model, const size_t size)
  {
    Model.SetXCellSizes().resize(size);
    Model.SetYCellSizes().resize(size);
    Model.SetZCellSizes().resize(size);

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

int main()
  {
    //for each block size we perform several runs and average the run time
    //to reduce the influence of other running programs
    const size_t nrunspersize = 5;
    std::string filename = "blocksize.time";
    //we create a calculator and implementation object manually
    //we do not use the factory function to have maximum control
    boost::shared_ptr<jif3D::ThreeDGravMagCalculator<jif3D::ThreeDGravityModel> > Calculator;
    boost::shared_ptr<jif3D::TensorCudaGravityImp> Implementation(
        new jif3D::TensorCudaGravityImp);

    //we test for a number of different blocksizes
    //64 and 256 are recommended values from the NVidia documentation
    std::vector<double> blocksizes;
    blocksizes += 64, 76, 90, 107, 128, 152, 181, 192, 215, 256;

    const size_t nruns = blocksizes.size();
    std::ofstream outfile(filename.c_str());
    std::cout << " Starting calculations. " << std::endl;
    //go through the different block sizes
    for (size_t i = 0; i < nruns; ++i)
      {
        //we use a fixed model size with modelsize cells in each spatial direction
        const size_t modelsize = 80;
        //to show that something is happening we print the current block size to the screent
        std::cout << "Blocksize: " << blocksizes.at(i) << std::endl;
        jif3D::ThreeDGravityModel GravityTest;
        //set the block size in the implementation object
        Implementation->SetCUDABlockSize(blocksizes.at(i));
        //and assemble the calculator object
        Calculator = boost::shared_ptr<jif3D::ThreeDGravMagCalculator<jif3D::ThreeDGravityModel> >(
            new jif3D::MinMemGravMagCalculator<jif3D::ThreeDGravityModel>(Implementation));
        double rawruntime = 0.0;
        //now we perform several runs and measure the time
        for (size_t j = 0; j < nrunspersize; ++j)
          {
            rawruntime = 0.0;
            //create a random test model
            MakeTestModel(GravityTest, modelsize);
            //save the time when we started
            boost::posix_time::ptime firststarttime =
                boost::posix_time::microsec_clock::local_time();
            //perform the calculation
            jif3D::rvec gravmeas(Calculator->Calculate(GravityTest));
            //and the time we stop
            boost::posix_time::ptime firstendtime =
                boost::posix_time::microsec_clock::local_time();
            //add to the time for averaging
            rawruntime += (firstendtime - firststarttime).total_microseconds();
          }
        //divide by number of runs and output the average
        rawruntime /= nrunspersize;
        outfile << blocksizes.at(i) << " " << rawruntime << std::endl;
      }
    //end of program
  }
