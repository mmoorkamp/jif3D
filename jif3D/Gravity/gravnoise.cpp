//============================================================================
// Name        : gravnoise.cpp
// Author      : Max Moorkamp
// Version     :
// Copyright   : 2008, MM
//============================================================================

/*! \file gravnoise.cpp
 * Add random gaussian noise to a netcdf file with scalar gravity measurements.
 */

#include "ThreeDGravityModel.h"
#include "ReadWriteGravityData.h"
#include "../Global/FileUtil.h"
#include <iostream>
#include <ctime>
#include <boost/random/lagged_fibonacci.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/variate_generator.hpp>

int main()
  {
    //create the generator object for the random number generator
    boost::lagged_fibonacci607 generator(
        static_cast<unsigned int> (std::time(0)));

    std::vector<double> PosX, PosY, PosZ, Data, Error;

    //read in the netcdf file with the data
    std::string datafilename = jif3D::AskFilename( "Data Filename: ");
    jif3D::ReadScalarGravityMeasurements(datafilename, Data, PosX, PosY, PosZ, Error);

    //get the relative noise level
    double noiselevel;
    std::cout << "Relative noise level: ";
    std::cin >> noiselevel;

    const size_t nmeas = Data.size();
    Error.resize(nmeas);
    //create a gaussian distribution for each datum and draw a sample from it
    for (size_t i = 0; i < nmeas; ++i)
      {
        boost::normal_distribution<> dist(Data.at(i), fabs(Data.at(i) * noiselevel));
        boost::variate_generator<boost::lagged_fibonacci607&, boost::normal_distribution<> >
            Sample(generator, dist);
        Data.at(i) = Sample();
        Error.at(i) = Data.at(i) * noiselevel;
      }
    //write noisy data to a file
    jif3D::SaveScalarGravityMeasurements(datafilename+".noise.nc", Data, PosX, PosY, PosZ, Error);
  }
