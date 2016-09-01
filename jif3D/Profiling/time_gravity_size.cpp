#ifdef HAVEHPX
#include <hpx/hpx.hpp>
#include <hpx/hpx_init.hpp>
#endif
#ifdef HAVEOPENMP
#include <omp.h>
#endif
#include <iostream>
#include <fstream>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/program_options.hpp>
#include "../Global/convert.h"
#include "../Gravity/ThreeDGravityModel.h"
#include "../GravMag/ThreeDGravMagCalculator.h"
#include "../GravMag/FullSensitivityGravMagCalculator.h"
#include "../GravMag/MinMemGravMagCalculator.h"
#include "../GravMag/DiskGravMagCalculator.h"
#include "../Gravity/ThreeDGravityFactory.h"

/*! \file time_gravity_size.cpp
 * This program can be used to benchmark the different methods for calculating
 * scalar and tensor gravity data as a function of model size. Different methods
 * can be set by options on the command line. To get a list of options type
 * time_gravity_size --help.
 */

void MakeTestModel(jif3D::ThreeDGravityModel &Model, const size_t size)
  {
    Model.SetMeshSize(size, size, size);

    for (size_t i = 0; i < size; ++i) // set the values of the inner cells
      {
        Model.SetXCellSizes()[i] = rand() % 10000 + 1000;
        Model.SetYCellSizes()[i] = rand() % 10000 + 1000;
        Model.SetZCellSizes()[i] = 500;
      }
    //fill the grid with some random values that are in a similar range
    // as densities
    std::generate_n(Model.SetDensities().origin(), Model.SetDensities().num_elements(),
        []() -> double
          { return 0.1 + double(rand() % 1000) / 300.0;});

    const size_t nmeas = 30;
    for (size_t i = 0; i < nmeas; ++i)
      Model.AddMeasurementPoint(rand() % 50000 + 2e4, rand() % 50000 + 2e4, 0.0);
    std::vector<double> bg_dens =
      { 1.0, 1.0, 5.0, 5.0 };
    std::vector<double> bg_thick =
      { 200.0, 300.0, 3500.0, 1000.0 };
    Model.SetBackgroundDensities(bg_dens);
    Model.SetBackgroundThicknesses(bg_thick);
  }

namespace po = boost::program_options;
int caching = 0;
po::options_description desc("Allowed options");

int hpx_main(boost::program_options::variables_map& vm)
  {

    if (vm.count("help"))
      {

        std::cout << desc << "\n";
#ifdef HAVEHPX
        return hpx::finalize();
#endif
        return 1;
      }
    const size_t nruns = 50;
    const size_t nrunspersize = 5;
    std::string filename;
    bool wantcuda = false;
    boost::shared_ptr<jif3D::ThreeDGravMagCalculator<jif3D::ThreeDGravityModel> > Calculator;

    if (vm.count("gpu"))
      {
        filename = "gpu_";
        std::cout << "Using GPU" << "\n";
        wantcuda = true;
      }
    else
      {
        filename = "cpu_";
#ifdef HAVEOPENMP
        if (vm.count("threads"))
          {
            omp_set_num_threads(vm["threads"].as<int>());
            filename += jif3D::stringify(vm["threads"].as<int>());
          }
        else
          {
            filename += jif3D::stringify(omp_get_max_threads());
          }
#endif
#ifdef HAVEHPX
        std::vector<hpx::naming::id_type> localities = hpx::find_all_localities();
        const size_t nthreads = hpx::get_num_worker_threads();
        const size_t nlocs = localities.size();
        filename += "l" + jif3D::stringify(nlocs) + "t" + jif3D::stringify(nthreads);
#endif
      }

    if (vm.count("ftg"))
      {
        filename += "ftg.time";
        switch (caching)
          {
        case 2:
          Calculator =
              jif3D::CreateGravityCalculator<
                  jif3D::FullSensitivityGravMagCalculator<jif3D::ThreeDGravityModel> >::MakeTensor(
                  wantcuda);
          break;
        case 1:
          Calculator = jif3D::CreateGravityCalculator<
              jif3D::DiskGravMagCalculator<jif3D::ThreeDGravityModel> >::MakeTensor(
              wantcuda);
          break;
        case 0:
          Calculator = jif3D::CreateGravityCalculator<
              jif3D::MinMemGravMagCalculator<jif3D::ThreeDGravityModel> >::MakeTensor(
              wantcuda);
          break;

          }
      }
    else
      {
        filename += "scalar.time";
        switch (caching)
          {
        case 2:
          Calculator =
              jif3D::CreateGravityCalculator<
                  jif3D::FullSensitivityGravMagCalculator<jif3D::ThreeDGravityModel> >::MakeScalar(
                  wantcuda);
          break;
        case 1:
          Calculator = jif3D::CreateGravityCalculator<
              jif3D::DiskGravMagCalculator<jif3D::ThreeDGravityModel> >::MakeScalar(
              wantcuda);
          break;
        case 0:
          Calculator = jif3D::CreateGravityCalculator<
              jif3D::MinMemGravMagCalculator<jif3D::ThreeDGravityModel> >::MakeScalar(
              wantcuda);
          break;

          }
      }

    std::ofstream outfile(filename.c_str());
    std::cout << " Starting calculations. " << std::endl;
    // we calculate gravity data for a number of different grid sizes
    for (size_t i = 0; i < nruns; ++i)
      {
        const size_t modelsize = (i + 1) * 2;
        std::cout << "Current model size: " << pow(modelsize, 3) << std::endl;
        jif3D::ThreeDGravityModel GravityTest;

        double rawruntime = 0.0;
        double cachedruntime = 0.0;
        //for each grid size we perform several runs and average the run time
        //to reduce the influence from other processes running on the system
        for (size_t j = 0; j < nrunspersize; ++j)
          {
            rawruntime = 0.0;
            cachedruntime = 0.0;
            MakeTestModel(GravityTest, modelsize);

            boost::posix_time::ptime firststarttime =
                boost::posix_time::microsec_clock::local_time();
            jif3D::rvec gravmeas(Calculator->Calculate(GravityTest));

            boost::posix_time::ptime firstendtime =
                boost::posix_time::microsec_clock::local_time();
            rawruntime += (firstendtime - firststarttime).total_microseconds();
            //if we want to compare to the time using caching we calculate
            //a cached result right after the original run
            if (caching > 0)
              {
                boost::posix_time::ptime secondstarttime =
                    boost::posix_time::microsec_clock::local_time();
                jif3D::rvec gravmeas2(Calculator->Calculate(GravityTest));
                boost::posix_time::ptime secondendtime =
                    boost::posix_time::microsec_clock::local_time();

                cachedruntime += (secondendtime - secondstarttime).total_microseconds();
              }
          }
        rawruntime /= nrunspersize;
        cachedruntime /= nrunspersize;
        outfile << modelsize * modelsize * modelsize << " " << rawruntime << " "
            << cachedruntime << std::endl;
      }
#ifdef HAVEHPX
    return hpx::finalize();
#endif
    return 0;
  }

int main(int argc, char* argv[])
  {

    desc.add_options()("help", "produce help message")("scalar",
        "Perform scalar calculation [default]")("ftg", "Perform FTG calculation ")("cpu",
        "Perform calculation on CPU [default]")("gpu", "Perform calculation on GPU")(
        "cachetype", po::value<int>(&caching)->default_value(0),
        "0 = no caching, 1 = disk, 2 = memory")("threads", po::value<int>(),
        "The number of openmp threads");
#ifdef HAVEHPX
    return hpx::init(desc, argc, argv);
#else
//set up the command line options
    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);
    return hpx_main(vm);
#endif
  }
