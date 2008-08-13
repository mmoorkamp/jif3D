//============================================================================
// Name        : gravnoise.cpp
// Author      : Max Moorkamp
// Version     :
// Copyright   : 2008, MM
//============================================================================

#include <iostream>
#include "ThreeDGravityModel.h"
#include <boost/random/lagged_fibonacci.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/variate_generator.hpp>
#include <ctime>

int main(int argc, char *argv[])
  {
    boost::lagged_fibonacci607 generator(
        static_cast<unsigned int> (std::time(0)));

    jiba::ThreeDGravityModel::tScalarMeasVec Data;
    jiba::ThreeDGravityModel::tMeasPosVec PosX, PosY, PosZ;

    std::string datafilename;
    std::cout << "Data Filename: ";
    std::cin >> datafilename;
    jiba::ReadScalarGravityMeasurements(datafilename, Data, PosX, PosY, PosZ);

    double noiselevel;
    std::cout << "Relative noise level: ";
    std::cin >> noiselevel;

    const size_t nmeas = Data.size();
    for (size_t i = 0; i < nmeas; ++i)
      {
        boost::normal_distribution<> dist(Data.at(i), Data.at(i) * noiselevel);
        boost::variate_generator<boost::lagged_fibonacci607&, boost::normal_distribution<> >
            Sample(generator, dist);
        Data.at(i) = Sample();
      }
    jiba::SaveScalarGravityMeasurements(datafilename+".noise.nc", Data, PosX, PosY, PosZ);
  }
