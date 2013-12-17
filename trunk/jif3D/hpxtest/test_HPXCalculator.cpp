//============================================================================
// Name        : test_ReadWriteX3D.cpp
// Author      : Feb 17, 2009
// Version     :
// Copyright   : 2009, mmoorkamp
//============================================================================

#include <hpx/config.hpp>
#include <hpx/hpx_init.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/program_options.hpp>
#include <boost/program_options/config.hpp>
#include "../Global/NumUtil.h"
#include "HPXMTCalculator.h"
#include "../MT/X3DModel.h"
#include "../MT/X3DMTCalculator.h"
#include "../MT/MTEquations.h"
#include "../MT/MT2DForward.h"
#include "../MT/ReadWriteImpedances.h"

namespace po = boost::program_options;

int hpx_main(int argc, char* argv[])
  {

    po::options_description desc("General options");
    desc.add_options()("help", "produce help message")("tempdir",
        po::value<std::string>(),
        "The name of the directory to store temporary files in");
    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    boost::filesystem::path TempDir = boost::filesystem::current_path();
    if (vm.count("tempdir"))
      {
        TempDir = vm["tempdir"].as<std::string>();
        if (!boost::filesystem::is_directory(TempDir))
          {
            std::cerr << TempDir.string() << " is not a directory or does not exist ! \n";
            return 500;
          }
      }

    //create a 3D version of a 1D model
    const size_t xsize = 20;
    const size_t ysize = 20;
    const size_t zsize = 10;
    const size_t nbglayers = 5;
    jif3D::X3DModel Model;

    Model.SetZCellSizes().resize(boost::extents[zsize]);

    Model.SetConductivities().resize(boost::extents[xsize][ysize][zsize]);
    std::vector<double> bg_thicknesses(nbglayers), bg_conductivities(nbglayers);

    const double deltax = 100.0;
    const double deltay = 100.0;
    const double deltaz = 100.0;
    const double freq = 1.0;
    const size_t nfreq = 32;
    const double cond = 0.01;
    Model.SetHorizontalCellSize(deltax, deltay, xsize, ysize);

    std::fill_n(Model.SetZCellSizes().origin(), zsize, deltaz);
    std::fill_n(Model.SetConductivities().origin(), xsize * ysize * zsize, cond);
    for (size_t i = 0; i < xsize * ysize * zsize; ++i)
      {
        *(Model.SetConductivities().origin() + i) *= i * 1.01;
      }
    std::fill_n(bg_conductivities.begin(), nbglayers, cond);
    bg_conductivities.back() *= 1.001;
    std::fill_n(bg_thicknesses.begin(), nbglayers, 200.0);

    Model.SetBackgroundConductivities(bg_conductivities);
    Model.SetBackgroundThicknesses(bg_thicknesses);
    for (size_t i = 0; i < nfreq; ++i)
      {
        Model.SetFrequencies().push_back(freq * i + 1);
      }
    jif3D::X3DMTCalculator X3DCalculator(TempDir);
    jif3D::HPXMTCalculator HPXCalculator(TempDir);

    for (size_t i = 1; i < xsize / 2; ++i)
      for (size_t j = 1; j < ysize / 2; ++j)
        {
          Model.AddMeasurementPoint(Model.GetXCoordinates()[i] + deltax / 2.0,
              Model.GetYCoordinates()[j] + deltay / 2.0, 0.0);
        }
    std::cout << "X3D " << std::endl;
    boost::posix_time::ptime x3dstarttime =
        boost::posix_time::microsec_clock::local_time();
    jif3D::rvec X3DImpedance = X3DCalculator.Calculate(Model);
    boost::posix_time::ptime x3dendtime = boost::posix_time::microsec_clock::local_time();
    std::cout << "HPX " << std::endl;
    boost::posix_time::ptime hpxstarttime =
        boost::posix_time::microsec_clock::local_time();
    jif3D::rvec HPXImpedance = HPXCalculator.Calculate(Model);
    boost::posix_time::ptime hpxendtime = boost::posix_time::microsec_clock::local_time();
    std::complex<double> HsImp = jif3D::ImpedanceHalfspace(freq, cond);
    const double prec = 0.05;
    const size_t nsites = HPXImpedance.size() / 8;
    for (size_t i = 0; i < nsites; ++i)
      {
        if (std::abs(
            (HPXImpedance(i * 8 + 2) - X3DImpedance(i * 8 + 2)) / X3DImpedance(i * 8 + 2))
            > 0.01
            || std::abs(
                (HPXImpedance(i * 8 + 3) - X3DImpedance(i * 8 + 3))
                    / X3DImpedance(i * 8 + 3)) > 0.01
            || std::abs(
                (HPXImpedance(i * 8 + 4) - X3DImpedance(i * 8 + 4))
                    / X3DImpedance(i * 8 + 4)) > 0.01
            || std::abs(
                (HPXImpedance(i * 8 + 5) - X3DImpedance(i * 8 + 5))
                    / X3DImpedance(i * 8 + 5)) > 0.01)
          {
            std::cout << HPXImpedance(i * 8 + 2) << " " << X3DImpedance(i * 8 + 2)
                << "     " << HPXImpedance(i * 8 + 3) << " " << X3DImpedance(i * 8 + 3)
                << std::endl;
            std::cout << HPXImpedance(i * 8 + 4) << " " << X3DImpedance(i * 8 + 4)
                << "     " << HPXImpedance(i * 8 + 5) << " " << X3DImpedance(i * 8 + 5)
                << std::endl;
          }
      }
    std::cout << " Times: " << std::endl;
    std::cout << "x3d: " << (x3dendtime - x3dstarttime).total_microseconds() << std::endl;
    std::cout << "hpx: " << (hpxendtime - hpxstarttime).total_microseconds() << std::endl;
    // Any HPX application logic goes here...
    return hpx::finalize();
    return 1;
  }

int main(int argc, char* argv[])
  {

    // Initialize HPX, run hpx_main as the first HPX thread, and
    // wait for hpx::finalize being called.

    return hpx::init(argc, argv);

  }

