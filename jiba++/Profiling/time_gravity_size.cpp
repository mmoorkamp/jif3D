#include <iostream>
#include <boost/date_time/posix_time/posix_time.hpp>
#include "../Gravity/ThreeDGravityModel.h"
#include "../Gravity/FullSensitivityGravityCalculator.h"
#include "../Gravity/WaveletCompressedGravityCalculator.h"

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
        for (size_t k = 0; k < size; ++k)
          {
            Model.SetDensities()[i][j][k] = double(rand() % 1000) / 300.0;
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
    const size_t nruns = 100;
    const size_t nrunspersize = 5;
    for (size_t i = 0; i < nruns; ++i)
      {
        const size_t modelsize = (i + 1) * 2;

        jiba::ThreeDGravityModel GravityTest;

        double rawruntime = 0.0;
        double cachedruntime = 0.0;
        for (size_t j = 0; j < nrunspersize; ++j)
          {
            rawruntime = 0.0;
            cachedruntime = 0.0;
            MakeRandomModel(GravityTest, modelsize);

            boost::shared_ptr<jiba::FullSensitivityGravityCalculator>
                Calculator = jiba::CreateGravityCalculator<
                    jiba::FullSensitivityGravityCalculator>::MakeScalar(false);
            boost::posix_time::ptime firststarttime =
                        boost::posix_time::microsec_clock::local_time();
            jiba::rvec gravmeas(Calculator->Calculate(GravityTest));

            boost::posix_time::ptime firstendtime =
                boost::posix_time::microsec_clock::local_time();

            boost::posix_time::ptime secondstarttime =
                boost::posix_time::microsec_clock::local_time();
            jiba::rvec gravmeas2(Calculator->Calculate(GravityTest));
            boost::posix_time::ptime secondendtime =
                boost::posix_time::microsec_clock::local_time();
            rawruntime += (firstendtime
                - firststarttime).total_microseconds();
            cachedruntime += (secondendtime
                - secondstarttime).total_microseconds();
          }
        rawruntime /= nrunspersize;
        cachedruntime /= nrunspersize;
        std::cout << modelsize * modelsize * modelsize << " " << rawruntime << " " << cachedruntime << std::endl;

      }

  }
