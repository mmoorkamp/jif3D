//============================================================================
// Name        : makegravmesh.cpp
// Author      : Max Moorkamp
// Version     :
// Copyright   : 2008, MM
//============================================================================

#ifdef HAVEHPX
#include <hpx/config.hpp>
#include <hpx/hpx_init.hpp>
#include <hpx/include/actions.hpp>
#include <hpx/include/util.hpp>
#include <hpx/include/lcos.hpp>
#include <hpx/include/iostreams.hpp>
#endif

#ifdef HAVEOPENMP
#include <omp.h>
#endif

#include <iostream>
#include <string>
#include <cstdlib>

#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/serialization/serialization.hpp>
#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>
#include "../Global/FileUtil.h"
#include "../Global/convert.h"
#include "../MT/X3DModel.h"
#include "../MT/X3DMTCalculator.h"

struct realinfo
  {
  double bg_conductivity;
  double phase1cond;
  double phase2cond;
  double phase1frac;
  std::string tempdir;
  std::string x3dname;

  //! Provide serialization to be able to store objects
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version)
    {
      ar & bg_conductivity;
      ar & phase1cond;
      ar & phase2cond;
      ar & phase1frac;
      ar & tempdir;
      ar & x3dname;
    }
  realinfo(double bgc, double p1c, double p2c, double p1f, std::string td,
      std::string x3d) :
      bg_conductivity(bgc), phase1cond(p1c), phase2cond(p2c), phase1frac(p1f), tempdir(
          td), x3dname(x3d)
    {

    }
  realinfo() :
      bg_conductivity(-1.0), phase1cond(-1.0), phase2cond(-1.0), phase1frac(-1.0), tempdir(), x3dname()
    {

    }
  };

jif3D::rvec CalcRealization(jif3D::X3DModel Model, realinfo Info);

#ifdef HAVEHPX
HPX_PLAIN_ACTION(CalcRealization, Calc_action);
#endif

jif3D::rvec CalcRealization(jif3D::X3DModel Model, realinfo Info)
  {
    const size_t nx = Model.GetXCoordinates().num_elements();
    const size_t ny = Model.GetYCoordinates().num_elements();
    const size_t nz = Model.GetZCoordinates().num_elements();

    for (size_t i = 0; i < nx; ++i)
      {
        for (size_t j = 0; j < ny; ++j)
          {
            Model.SetConductivities()[i][j][0] = Info.bg_conductivity;
            for (size_t k = 1; k < nz; ++k)
              {
                if (drand48() < Info.phase1frac)
                  {
                    Model.SetConductivities()[i][j][k] = Info.phase1cond;
                  }
                else
                  {
                    Model.SetConductivities()[i][j][k] = Info.phase2cond;
                  }
              }
          }
      }

    jif3D::X3DMTCalculator Calculator(Info.tempdir, Info.x3dname);
    jif3D::rvec Impedances(Calculator.Calculate(Model));
    return Impedances;
  }

namespace po = boost::program_options;

int hpx_main(po::variables_map& vm)
  {
#ifdef HAVEHPX
    using hpx::cout;
#else
    using std::cout;
#endif

    using std::cin;

#ifdef HAVEOPENMP
    if (vm.count("threads"))
      {
        omp_set_num_threads(vm["threads"].as<int>());
      }
#endif

    double topthick = -1.0;
    std::string x3dname("x3d");
    std::string tempdir;

    if (vm.count("topthick"))
      {
        topthick = vm["topthick"].as<double>();
      }

    if (vm.count("tempdir"))
      {
        tempdir = vm["tempdir"].as<std::string>();
      }
    else
      {
        tempdir = boost::filesystem::current_path().native();
      }

    if (vm.count("x3dname"))
      {
        x3dname = vm["x3dname"].as<std::string>();
      }
    jif3D::X3DModel Model;
    size_t nx, ny, nz;
    double deltax, deltay, deltaz;
    //first find out the basic mesh parameters
    //the number of cells in each coordinate direction
    cout << "Number of cells Nx: ";
    cin >> nx;
    cout << "Number of cells Ny: ";
    cin >> ny;
    cout << "Number of cells Nz: ";
    cin >> nz;
    //and the size of the cells in each direction
    cout << "Cell size x [m]: ";
    cin >> deltax;
    cout << "Cell size y [m]: ";
    cin >> deltay;
    cout << "Cell size z [m]: ";
    cin >> deltaz;
    //set the cell sizes and allocate memory for the mesh
    Model.SetMeshSize(nx, ny, nz);
    Model.SetHorizontalCellSize(deltax, deltay, nx, ny);
    std::fill_n(Model.SetZCellSizes().begin(), nz, deltaz);
    if (topthick > 0.0)
      {
        Model.SetZCellSizes()[0] = topthick;
      }
    //ask for a conductivity to fill the mesh with
    double bg_conductivity = 1.0;
    cout << "Background Conductivity [S/m] : ";
    cin >> bg_conductivity;
    double phase1cond, phase2cond, phase1frac;
    cout << "Conductivity Phase 1 [S/m] : ";
    cin >> phase1cond;
    cout << "Conductivity Phase 2 [S/m] : ";
    cin >> phase2cond;
    cout << "Fraction Phase 1: ";
    cin >> phase1frac;
    srand48(time(0));
    double frequency;
    cout << "Frequency [Hz] : ";
    cin >> frequency;
    double posx, posy, posz = 0;
    posx = (deltax * nx) / 2.0;
    posy = (deltay * ny) / 2.0;
    Model.SetFrequencies().assign(1, frequency);
    Model.AddMeasurementPoint(posx, posy, posz);

    //ask for a filename to write the mesh to
    std::string OutFilename = jif3D::AskFilename("Outfile name: ", false);

    //fill the background
    std::vector<double> bg_thicknesses(Model.GetZCellSizes().size(), deltaz);
    std::vector<double> bg_conductivities(Model.GetZCellSizes().size(), bg_conductivity);
    Model.SetBackgroundThicknesses(bg_thicknesses);
    Model.SetBackgroundConductivities(bg_conductivities);
    size_t nrealmax;
    cout << "Realizations: ";
    cin >> nrealmax;

    boost::posix_time::ptime starttime = boost::posix_time::microsec_clock::local_time();

    std::ofstream zxy((OutFilename + "_zxy.out").c_str());
    std::ofstream zyx((OutFilename + "_zyx.out").c_str());
    std::ofstream zxx((OutFilename + "_zxx.out").c_str());
    std::ofstream zyy((OutFilename + "_zyy.out").c_str());

#ifdef HAVEHPX
    using hpx::lcos::unique_future;
    using hpx::async;
    using hpx::wait_all;
    std::vector<unique_future<jif3D::rvec> > ImplResult;
    ImplResult.reserve(nrealmax);
    Calc_action CalcImpl;
    std::vector<hpx::naming::id_type> localities = hpx::find_all_localities();

    cout << "Found " << localities.size() << " localities\n ";

    for (size_t nreal = 0; nreal < nrealmax; ++nreal)
      {
        hpx::naming::id_type const locality_id = localities.at(nreal % localities.size());
        ImplResult.push_back(
            async(CalcImpl, locality_id, Model, realinfo(bg_conductivity, phase1cond, phase2cond,
                    phase1frac, tempdir, x3dname)));
      }
    wait_all(ImplResult);

    for (size_t nreal = 0; nreal < nrealmax; ++nreal)
      {
        jif3D::rvec Impedances = ImplResult[nreal].get();
        std::string realstring(jif3D::stringify(nreal));
        zxx << nreal << " " << Impedances(0) << " " << Impedances(1) << std::endl;
        zxy << nreal << " " << Impedances(2) << " " << Impedances(3) << std::endl;
        zyx << nreal << " " << Impedances(4) << " " << Impedances(5) << std::endl;
        zyy << nreal << " " << Impedances(6) << " " << Impedances(7) << std::endl;

      }
    boost::posix_time::ptime endtime = boost::posix_time::microsec_clock::local_time();
    double runtime = (endtime - starttime).total_seconds();
    cout << "Runtime: " << runtime << " s" << hpx::endl;
    return hpx::finalize();
#else
#pragma omp parallel for shared(Model, zxy, zyx, zxx, zyy)
    for (size_t nreal = 0; nreal < nrealmax; ++nreal)
      {

        jif3D::rvec Impedances = CalcRealization(Model,
            realinfo(bg_conductivity, phase1cond, phase2cond, phase1frac, tempdir,
                x3dname));
        jif3D::rvec Errors(Impedances.size(), 0.0);
        std::string realstring(jif3D::stringify(nreal));
#pragma omp critical(write_files)
          {
            //RealModel.WriteNetCDF(OutFilename + realstring + ".nc");
            //RealModel.WriteVTK(OutFilename + realstring + ".vtk");
            zxx << nreal << " " << Impedances(0) << " " << Impedances(1) << std::endl;
            zxy << nreal << " " << Impedances(2) << " " << Impedances(3) << std::endl;
            zyx << nreal << " " << Impedances(4) << " " << Impedances(5) << std::endl;
            zyy << nreal << " " << Impedances(6) << " " << Impedances(7) << std::endl;
          }
      }

    boost::posix_time::ptime endtime = boost::posix_time::microsec_clock::local_time();
    double runtime = (endtime - starttime).total_seconds();
    cout << "Runtime: " << runtime << " s\n";
    return 0;
#endif
  }

int main(int argc, char* argv[])
  {
    po::options_description desc("General options");
    desc.add_options()("help", "produce help message")("threads", po::value<int>(),
        "The number of openmp threads")("topthick", po::value<double>(),
        "Thickness of the top layer in m, if <= 0.0 the thickness will be the same as the z dimension of the other grid cells.")(
        "x3dname", po::value<std::string>(), "The name of the executable for x3d")(
        "tempdir", po::value<std::string>(),
        "The name of the directory where we store files for forward calculation");
#ifdef HAVEHPX
    return hpx::init(desc, argc, argv);
#else
    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);

    if (vm.count("help"))
      {
        std::cout << desc << "\n";
        return 1;
      }
    return hpx_main(vm);
#endif
  }
