//============================================================================
// Name        : gravnoise.cpp
// Author      : Max Moorkamp
// Version     :
// Copyright   : 2008, MM
//============================================================================

/*! \file gravnoise.cpp
 * Add random gaussian noise to a netcdf file with scalar gravity measurements.
 */

#include <iostream>
#include <ctime>
#include <boost/random/lagged_fibonacci.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/variate_generator.hpp>
#include "../Global/FileUtil.h"
#include "ThreeDGravityModel.h"
#include "ReadWriteGravityData.h"


int main()
  {
    //create the generator object for the random number generator
    boost::lagged_fibonacci607 generator(
        static_cast<unsigned int> (std::time(0)));

    jif3D::rvec Data;
    jif3D::ThreeDGravityModel::tMeasPosVec PosX, PosY, PosZ;

    //read in the netcdf file with the data
    std::string datafilename = jif3D::AskFilename( "Data Filename: ");
    jif3D::ReadScalarGravityMeasurements(datafilename, Data, PosX, PosY, PosZ);

    //get the relative noise level
    double noiselevel;
    std::cout << "Relative noise level: ";
    std::cin >> noiselevel;

    const size_t nmeas = Data.size();
    //create a gaussian distribution for each datum and draw a sample from it
    for (size_t i = 0; i < nmeas; ++i)
      {
        boost::normal_distribution<> dist(Data(i), fabs(Data(i) * noiselevel));
        boost::variate_generator<boost::lagged_fibonacci607&, boost::normal_distribution<> >
            Sample(generator, dist);
        Data(i) = Sample();
      }
    //write noisy data to a file
    jif3D::SaveScalarGravityMeasurements(datafilename+".noise.nc", Data, PosX, PosY, PosZ);
  }
