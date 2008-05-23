#include <iostream>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <ThreeDGravityModel.h>

int main()
  {
    boost::posix_time::ptime starttime =
        boost::posix_time::microsec_clock::local_time();

    jiba::ThreeDGravityModel GravityTest(true, true);
    //create a model with a random number of cells in each direction
    const size_t nhorcells = 100;
    const size_t nzcells = 20;
    GravityTest.SetXCellSizes().resize(boost::extents[nhorcells]);
    GravityTest.SetYCellSizes().resize(boost::extents[nhorcells]);
    GravityTest.SetZCellSizes().resize(boost::extents[nzcells]);
    for (size_t i = 0; i < nhorcells; ++i) // set the values of the inner cells
      {
        GravityTest.SetXCellSizes()[i] = rand() % 10000 + 1000;
        GravityTest.SetYCellSizes()[i] = rand() % 10000 + 1000;
      }
    for (size_t i = 0; i < nzcells; ++i)
      {
        GravityTest.SetZCellSizes()[i] = 500;
      }

    GravityTest.SetDensities().resize(boost::extents[nhorcells][nhorcells][nzcells]);
    for (size_t i = 0; i < nhorcells; ++i)
      for (size_t j = 0; j < nhorcells; ++j)
        GravityTest.SetDensities()[i][j][0] = 1.0;
    for (size_t i = 0; i < nhorcells; ++i)
      for (size_t j = 0; j < nhorcells; ++j)
        for (size_t k = 1; k < nzcells; ++k)
          GravityTest.SetDensities()[i][j][k] = 5.0;

    const size_t nmeas = 30;
    for (size_t i = 0; i < nmeas; ++i)
      GravityTest.AddMeasurementPoint(rand() % 50000 + 2e4, rand() % 50000
          + 2e4, 0.0);
    std::vector<double> bg_dens, bg_thick;
    bg_dens.push_back(1.0);
    bg_dens.push_back(1.0);
    bg_dens.push_back(5.0);
    bg_dens.push_back(5.0);
    bg_thick.push_back(200.0);
    bg_thick.push_back(300.0);
    bg_thick.push_back(3500.0);
    bg_thick.push_back(1000.0);
    GravityTest.SetBackgroundDensities(bg_dens);
    GravityTest.SetBackgroundThicknesses(bg_thick);
    jiba::ThreeDGravityModel::tScalarMeasVec
        gravmeas(GravityTest.CalcGravity());

    boost::posix_time::ptime endtime =
        boost::posix_time::microsec_clock::local_time();
    std::cout << "First Calculation Runtime (initial calculation) : "
        << endtime - starttime << std::endl;
    starttime = boost::posix_time::microsec_clock::local_time();
    jiba::ThreeDGravityModel::tScalarMeasVec
        gravmeas2(GravityTest.CalcGravity());
    endtime = boost::posix_time::microsec_clock::local_time();
    std::cout << "Second Calculation Runtime (stored sensitvities) : "
        << endtime - starttime << std::endl;
    for (size_t i = 0; i < nhorcells; ++i) // set the values of the inner cells
      {
        GravityTest.SetXCellSizes()[i] = rand() % 10000 + 1000;
        GravityTest.SetYCellSizes()[i] = rand() % 10000 + 1000;
      }
    starttime = boost::posix_time::microsec_clock::local_time();
    jiba::ThreeDGravityModel::tScalarMeasVec
        gravmeas3(GravityTest.CalcGravity());
    endtime = boost::posix_time::microsec_clock::local_time();
    std::cout << "Third Calculation Runtime (after geometry change) : "
        << endtime - starttime << std::endl;
    for (int i = 0; i < nmeas; ++i)
      std::cout << gravmeas.at(i) << " " << gravmeas2.at(i) << " " << gravmeas3.at(i) << std::endl; 

  }
