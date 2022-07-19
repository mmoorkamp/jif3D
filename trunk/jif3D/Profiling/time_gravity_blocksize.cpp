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
#include <boost/make_shared.hpp>

using namespace boost::assign;

/*! \file time_gravity_blocksize.cpp
 * This program shows the execution time for the gavity forward calculation
 * using CUDA for a fixed model size, but different CUDA thread block sizes.
 */

void MakeTestModel(jif3D::ThreeDGravityModel &Model, jif3D::GeneralData &Data,
    const size_t size)
  {
    Model.SetMeshSize(size, size, size);
    jif3D::ThreeDModelBase::t3DModelDim XCS(size), YCS(size), ZCS(size, 500);
    Model.SetZCellSizes(ZCS);
    std::generate(XCS.begin(), XCS.end(), []() -> double
      { return rand() % 10000 + 1000;});
    std::generate(YCS.begin(), YCS.end(), []() -> double
      { return rand() % 10000 + 1000;});
    Model.SetXCellSizes(XCS);
    Model.SetYCellSizes(YCS);
    //fill the grid with some random values that are in a similar range
    // as densities
    std::generate_n(Model.SetDensities().origin(), Model.SetDensities().num_elements(),
        []() -> double
          { return 0.1 + double(rand() % 1000) / 300.0;});

    const size_t nmeas = 30;
    for (size_t i = 0; i < nmeas; ++i)
      Data.AddMeasurementPoint(rand() % 50000 + 2e4, rand() % 50000 + 2e4, 0.0);
    std::vector<double> bg_dens =
      { 1.0, 1.0, 5.0, 5.0 };
    std::vector<double> bg_thick =
      { 200.0, 300.0, 3500.0, 1000.0 };
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
    boost::shared_ptr<jif3D::ThreeDGravMagCalculator<jif3D::TensorGravityData> > Calculator;
    boost::shared_ptr<jif3D::TensorCudaGravityImp> Implementation(
        new jif3D::TensorCudaGravityImp);

    //we test for a number of different blocksizes
    //64 and 256 are recommended values from the NVidia documentation
    std::vector<double> blocksizes;
    blocksizes += 64, 76, 90, 107, 128, 152, 181, 192, 215, 256;

    typename jif3D::TensorGravityData Data;


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
        Calculator =
            boost::make_shared<jif3D::MinMemGravMagCalculator<jif3D::TensorGravityData> >(
                Implementation);
        double rawruntime = 0.0;
        //now we perform several runs and measure the time
        for (size_t j = 0; j < nrunspersize; ++j)
          {
            rawruntime = 0.0;
            //create a random test model
            MakeTestModel(GravityTest, Data, modelsize);
            //save the time when we started
            boost::posix_time::ptime firststarttime =
                boost::posix_time::microsec_clock::local_time();
            //perform the calculation
            jif3D::rvec gravmeas(Calculator->Calculate(GravityTest, Data));
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
