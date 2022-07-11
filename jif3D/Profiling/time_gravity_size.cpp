#ifdef HAVEHPX
#include <hpx/hpx.hpp>
#include <hpx/hpx_init.hpp>
#include <hpx/hpx_init_params.hpp>
#include <hpx/runtime_distributed.hpp>
#include <hpx/modules/program_options.hpp>
namespace po = hpx::program_options;

#endif
#ifdef HAVEOPENMP
#include <omp.h>
#include <boost/program_options.hpp>
namespace po = boost::program_options;

#endif
#include <iostream>
#include <fstream>
#include <boost/date_time/posix_time/posix_time.hpp>
#include "../Global/convert.h"
#include "../Gravity/ThreeDGravityModel.h"
#include "../GravMag/ThreeDGravMagCalculator.h"
#include "../GravMag/FullSensitivityGravMagCalculator.h"
#include "../GravMag/MinMemGravMagCalculator.h"
#include "../GravMag/DiskGravMagCalculator.h"
#include "../Gravity/ThreeDGravityFactory.h"
#include "../Gravity/ScalarGravityData.h"
#include "../Gravity/TensorGravityData.h"

/*! \file time_gravity_size.cpp
 * This program can be used to benchmark the different methods for calculating
 * scalar and tensor gravity data as a function of model size. Different methods
 * can be set by options on the command line. To get a list of options type
 * time_gravity_size --help.
 */

int caching = 0;
po::options_description desc("Allowed options");

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

template<class CalculatorType>
void RunCalculation(CalculatorType &Calculator, const std::string &filename)
  {
    const size_t nruns = 50;
    const size_t nrunspersize = 5;
    std::ofstream outfile(filename.c_str());
    typename CalculatorType::DataType Data;
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
            MakeTestModel(GravityTest, Data, modelsize);

            boost::posix_time::ptime firststarttime =
                boost::posix_time::microsec_clock::local_time();
            jif3D::rvec gravmeas(Calculator.Calculate(GravityTest, Data));

            boost::posix_time::ptime firstendtime =
                boost::posix_time::microsec_clock::local_time();
            rawruntime += (firstendtime - firststarttime).total_microseconds();
            //if we want to compare to the time using caching we calculate
            //a cached result right after the original run
            if (caching > 0)
              {
                boost::posix_time::ptime secondstarttime =
                    boost::posix_time::microsec_clock::local_time();
                jif3D::rvec gravmeas2(Calculator.Calculate(GravityTest, Data));
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
  }

int hpx_main(po::variables_map &vm)
  {

    if (vm.count("help"))
      {

        std::cout << desc << "\n";
#ifdef HAVEHPX
        return hpx::finalize();
#endif
        return 1;
      }
    std::string filename;
    bool wantcuda = false;

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
        boost::shared_ptr<jif3D::ThreeDGravMagCalculator<jif3D::TensorGravityData> > Calculator;
        filename += "ftg.time";
        switch (caching)
          {
        case 2:
          Calculator =
              jif3D::CreateGravityCalculator<
                  jif3D::FullSensitivityGravMagCalculator<jif3D::TensorGravityData> >::MakeTensor(
                  wantcuda);
          break;
        case 1:
          Calculator = jif3D::CreateGravityCalculator<
              jif3D::DiskGravMagCalculator<jif3D::TensorGravityData> >::MakeTensor(
              wantcuda);
          break;
        case 0:
          Calculator = jif3D::CreateGravityCalculator<
              jif3D::MinMemGravMagCalculator<jif3D::TensorGravityData> >::MakeTensor(
              wantcuda);
          break;

          }
        RunCalculation(*Calculator, filename);
      }
    else
      {
        boost::shared_ptr<jif3D::ThreeDGravMagCalculator<jif3D::ScalarGravityData> > Calculator;
        filename += "scalar.time";
        switch (caching)
          {
        case 2:
          Calculator =
              jif3D::CreateGravityCalculator<
                  jif3D::FullSensitivityGravMagCalculator<jif3D::ScalarGravityData> >::MakeScalar(
                  wantcuda);
          break;
        case 1:
          Calculator = jif3D::CreateGravityCalculator<
              jif3D::DiskGravMagCalculator<jif3D::ScalarGravityData> >::MakeScalar(
              wantcuda);
          break;
        case 0:
          Calculator = jif3D::CreateGravityCalculator<
              jif3D::MinMemGravMagCalculator<jif3D::ScalarGravityData> >::MakeScalar(
              wantcuda);
          break;

          }
        RunCalculation(*Calculator, filename);
      }

#ifdef HAVEHPX
    return hpx::finalize();
#endif
    return 0;
  }

int main(int argc, char *argv[])
  {

    desc.add_options()("help", "produce help message")("scalar",
        "Perform scalar calculation [default]")("ftg", "Perform FTG calculation ")("cpu",
        "Perform calculation on CPU [default]")("gpu", "Perform calculation on GPU")(
        "cachetype", po::value<int>(&caching)->default_value(0),
        "0 = no caching, 1 = disk, 2 = memory")("threads", po::value<int>(),
        "The number of openmp threads");
#ifdef HAVEHPX
    hpx::init_params initparms;
    initparms.desc_cmdline = desc;
    return hpx::init(argc, argv, initparms);
#endif

//set up the command line options
    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);
    return hpx_main(vm);

  }
