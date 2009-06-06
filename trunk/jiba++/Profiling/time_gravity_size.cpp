#include <omp.h>
#include <iostream>
#include <fstream>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/program_options.hpp>
#include "../Global/convert.h"
#include "../Gravity/ThreeDGravityModel.h"
#include "../Gravity/ThreeDGravityCalculator.h"
#include "../Gravity/FullSensitivityGravityCalculator.h"
#include "../Gravity/MinMemGravityCalculator.h"

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

namespace po = boost::program_options;
int main(int ac, char* av[])
  {

    po::options_description desc("Allowed options");
    desc.add_options()("help", "produce help message")("scalar",
        "Perform scalar calculation [default]")("ftg",
        "Perform FTG calculation ")("cpu",
        "Perform calculation on CPU [default]")("gpu",
        "Perform calculation on GPU")("cached", "Also do cached calculation")(
        "threads", po::value<int>(), "The number of openmp threads");

    po::variables_map vm;
    po::store(po::parse_command_line(ac, av, desc), vm);
    po::notify(vm);

    if (vm.count("help"))
      {
        std::cout << desc << "\n";
        return 1;
      }

    const size_t nruns = 50;
    const size_t nrunspersize = 5;
    std::string filename;
    bool wantcuda = false;
    bool wantcached = false;
    boost::shared_ptr<jiba::ThreeDGravityCalculator> Calculator;


    if (vm.count("gpu"))
      {
        filename = "gpu_";
        std::cout << "Using GPU" << "\n";
        wantcuda = true;
      }
    else
      {
        filename = "cpu_";
        if (vm.count("threads"))
          {
            omp_set_num_threads(vm["threads"].as<int>());
            filename += jiba::stringify(vm["threads"].as<int>());
          }
        else
          {
            filename += jiba::stringify(omp_get_max_threads());
          }
      }


    if (vm.count("cached"))
      {
        wantcached = true;
        filename += "cached_";
      }

    if (vm.count("ftg"))
      {
        filename += "ftg.time";
        if (wantcached)
          {
            Calculator = jiba::CreateGravityCalculator<
                jiba::FullSensitivityGravityCalculator>::MakeTensor(wantcuda);
          }
        else
          {
            Calculator = jiba::CreateGravityCalculator<
                jiba::MinMemGravityCalculator>::MakeTensor(wantcuda);

          }

      }
    else
      {
        filename += "scalar.time";
        if (wantcached)
          {
            Calculator = jiba::CreateGravityCalculator<
                jiba::FullSensitivityGravityCalculator>::MakeScalar(wantcuda);
          }
        else
          {
            Calculator = jiba::CreateGravityCalculator<
                jiba::MinMemGravityCalculator>::MakeScalar(wantcuda);
          }
      }

    std::ofstream outfile(filename.c_str());
    std::cout << " Starting calculations. " << std::endl;
    for (size_t i = 0; i < nruns; ++i)
      {
        const size_t modelsize = (i + 1) * 2;
        std::cout << "Current model size: " << pow(modelsize, 3) << std::endl;
        jiba::ThreeDGravityModel GravityTest;

        double rawruntime = 0.0;
        double cachedruntime = 0.0;
        for (size_t j = 0; j < nrunspersize; ++j)
          {
            rawruntime = 0.0;
            cachedruntime = 0.0;
            MakeTestModel(GravityTest, modelsize);

            boost::posix_time::ptime firststarttime =
                boost::posix_time::microsec_clock::local_time();
            jiba::rvec gravmeas(Calculator->Calculate(GravityTest));

            boost::posix_time::ptime firstendtime =
                boost::posix_time::microsec_clock::local_time();
            rawruntime += (firstendtime - firststarttime).total_microseconds();
            if (wantcached)
              {
                boost::posix_time::ptime secondstarttime =
                    boost::posix_time::microsec_clock::local_time();
                jiba::rvec gravmeas2(Calculator->Calculate(GravityTest));
                boost::posix_time::ptime secondendtime =
                    boost::posix_time::microsec_clock::local_time();

                cachedruntime
                    += (secondendtime - secondstarttime).total_microseconds();
              }
          }
        rawruntime /= nrunspersize;
        cachedruntime /= nrunspersize;
        outfile << modelsize * modelsize * modelsize << " " << rawruntime
            << " " << cachedruntime << std::endl;
      }

  }
