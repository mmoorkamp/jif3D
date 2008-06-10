#include <iostream>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <ThreeDGravityModel.h>

void MakeRandomModel(jiba::ThreeDGravityModel &Model, const size_t size)
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
        for (int k = 0; k < size; ++k)
          {
            Model.SetDensities()[i][j][k] = double(rand() % 1000)/300.0;
          }
    const size_t nmeas = 30;
    for (size_t i = 0; i < nmeas; ++i)
      Model.AddMeasurementPoint(rand() % 50000 + 2e4, rand() % 50000
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
    Model.SetBackgroundDensities(bg_dens);
    Model.SetBackgroundThicknesses(bg_thick);
  }

int main()
  {
    const size_t nruns = 100;
    for (size_t i = 0; i < nruns; ++i)
      {
        const size_t modelsize = (i+1) * 2;
        boost::posix_time::ptime firststarttime =
            boost::posix_time::microsec_clock::local_time();
        jiba::ThreeDGravityModel GravityTest(true, true);

        MakeRandomModel(GravityTest, modelsize);

        jiba::ThreeDGravityModel::tScalarMeasVec
            gravmeas(GravityTest.CalcGravity());

        boost::posix_time::ptime firstendtime =
            boost::posix_time::microsec_clock::local_time();
        
        boost::posix_time::ptime secondstarttime =
                   boost::posix_time::microsec_clock::local_time();
        jiba::ThreeDGravityModel::tScalarMeasVec
            gravmeas2(GravityTest.CalcGravity());
        boost::posix_time::ptime secondendtime =
                   boost::posix_time::microsec_clock::local_time();
        
        std::cout << modelsize * modelsize * modelsize << " "  
            << (firstendtime - firststarttime).total_microseconds() << " " << (secondendtime - secondstarttime).total_microseconds() << std::endl;

      }

  }
