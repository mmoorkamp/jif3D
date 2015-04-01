#ifdef HAVEHPX
#include <hpx/config.hpp>
#include <hpx/hpx_init.hpp>
//#include <hpx/include/iostreams.hpp>
#endif
#ifdef HAVEOPENMP
#include <omp.h>
#endif
#include <iostream>
#include <fstream>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/program_options.hpp>
#include "../Global/convert.h"
#include "../Tomo/TomographyCalculator.h"

void MakeTestModel(jif3D::ThreeDSeismicModel &Model, const size_t size)
  {
    const double cellsize = 100.0;
    Model.SetCellSize(cellsize, size, size, size);
    const double firstdepth = Model.GetZCoordinates()[0];
    const double bottomdepth = Model.GetZCoordinates()[size - 1];
    const double topvel = 1000.0;
    const double bottomvel = 6000.0;

    Model.SetSlownesses().resize(boost::extents[size][size][size]);
    for (size_t i = 0; i < size; ++i)
      for (size_t j = 0; j < size; ++j)
        for (size_t k = 0; k < size; ++k)
          {
            double Depth = Model.GetZCoordinates()[k];
            double Velocity = topvel
                + (Depth - firstdepth) * (bottomvel - topvel)
                    / (bottomdepth - firstdepth);
            Model.SetSlownesses()[i][j][k] = 1.0 / Velocity;
          }

    const double minx = cellsize * 1.5;
    const double miny = minx;
    const double maxx = cellsize * size * 0.9;
    const double maxy = maxx;
    const double measz = 50.0;
    const double sourcez = 650;
    const size_t nmeasx = 10;
    const size_t nmeasy = 10;
    const double deltax = (maxx - minx) / nmeasx;
    const double deltay = (maxy - miny) / nmeasy;
    for (size_t i = 0; i <= nmeasx; ++i)
      {
        for (size_t j = 0; j <= nmeasy; ++j)
          {
            Model.AddSource(minx + i * deltax, miny + j * deltay, sourcez);
            Model.AddMeasurementPoint(minx + i * deltax, miny + j * deltay, measz);
          }
      }

    const size_t nsource = Model.GetSourcePosX().size();
    const size_t nmeas = Model.GetMeasPosX().size();
    for (size_t i = 0; i < nmeas; ++i)
      {
        for (size_t j = 0; j < nsource; ++j)
          {
            if (j != i)

              {
                Model.AddMeasurementConfiguration(j, i);
              }
          }
      }

  }

namespace po = boost::program_options;
int caching = 0;

int hpx_main(boost::program_options::variables_map& vm)
  {

    const size_t nruns = 50;
    const size_t nrunspersize = 5;
    std::string filename;
    bool wantcuda = false;
    jif3D::TomographyCalculator Calculator;

#ifdef HAVEHPX
    filename = "tomo_hpx_";
#endif

#ifdef HAVEOPENMP
    filename = "tomo_openmp_";
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

    std::ofstream outfile(filename.c_str());
    std::cout << " Starting calculations. " << std::endl;
    // we calculate gravity data for a number of different grid sizes
    for (size_t i = 0; i < nruns; ++i)
      {
        const size_t modelsize = 10 + (i + 1) * 2;
        std::cout << "Current model size: " << pow(modelsize, 3) << std::endl;
        jif3D::ThreeDSeismicModel SeisTest;

        double rawruntime = 0.0;

        //for each grid size we perform several runs and average the run time
        //to reduce the influence from other processes running on the system
        for (size_t j = 0; j < nrunspersize; ++j)
          {

            MakeTestModel(SeisTest, modelsize);

            boost::posix_time::ptime firststarttime =
                boost::posix_time::microsec_clock::local_time();
            jif3D::rvec gravmeas(Calculator.Calculate(SeisTest));

            boost::posix_time::ptime firstendtime =
                boost::posix_time::microsec_clock::local_time();
            rawruntime += (firstendtime - firststarttime).total_microseconds();

          }
        rawruntime /= nrunspersize;

        outfile << modelsize * modelsize * modelsize << " " << rawruntime << std::endl;
      }
#ifdef HAVEHPX
    return hpx::finalize();
#endif
    return 0;
  }

int main(int argc, char* argv[])
  {
    po::options_description desc("Allowed options");
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
